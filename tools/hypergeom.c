/** @file functions.c Documented special functions module.
 *
 * Nils Schoenberg, 12.10.2017
 *
 * This module computes complex gamma functions, the gaussian hypergeometric function (Hypergeometric2F1) and various derived functions (bessel integrals)
 */
#include "common.h"
#include "hypergeom.h"
#include <math.h>

#include <time.h>
//Mathematical constants used in this file

 #define MATH_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679 /** < The first 100 digits of pi*/


//Variations of PI
 #define MATH_PI_3_2 pow(MATH_PI,1.5) /** < PI^(3/2) */
 // 5.5683279968317078452848179821188357020136243902832 
 #define MATH_2_PI 2*MATH_PI/** < PI^(3/2) */
 #define MATH_PI_HALF 0.5*MATH_PI/** < PI/2 */
 #define SQRT_2_PI sqrt(2*MATH_PI)/** < sqrt(2*PI) */
//2.5066282746310005024157652848110452530069867406099;


//Logarithms of PI
 #define ln2PI log(2*MATH_PI)/** < ln(2*PI) */
 #define ln2PI_05 0.5*log(2*MATH_PI)/** < ln(2*PI)/2 */
 //0.39908993417905752478250359150769595020993410292128
 #define lnPI_HALF_05 0.5*log(MATH_PI*0.5)/** < ln(PI/2)/2 */
 
//Logarithms
 #define LOG_2 log(2)/** < log(2) */
  
 
 // Constants definition for the Lanczos Approximation
const int NUM_GAMMA_COEFF = 8;
//double GAMMA_COEFF[NUM_GAMMA_COEFF] = { 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
const double GAMMA_COEFF[] = { 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };

 
 // Precision parameters
double BESSEL_INTEGRAL_MAX_L = 80.0;
double BESSEL_INTEGRAL_T_SWITCH = 0.95;//0.95;
double BESSEL_INTEGRAL_T_SWITCH_SAFE = 0.9;//0.9;
double TAYLOR_HYPERGEOM_ACCURACY = 1e-10; 
double HYPERGEOM_EPSILON = 1e-20;
const int HYPERGEOM_MAX_ITER = 100000;
double ABS_LIMIT = 1e100; /** < Maximum acceptable absolute magnitude for a number, whose square should (easily) fit double precision */
double INV_ABS_LIMIT = 1e-100;
//Currently recursion relations are avoided at all costs
 double A_B_MAX = 10000; /** < Maximum value of parameters a and b for which we don't use recursion relations */
 double B_C_MAX = 10000; /** < Maximum value of parameters b and c for which we don't use recursion relations */
 double A_B_IMAG_LIMIT = 40.0; /** < Maximum value of parameters a and b for which no overflow occurs*/
 double NU_IMAG_LIMIT = 150; /** < Maximum value of nu_imag such that cosh(nu_imag) or gamma(nu_imag) doesn't overflow */
 
 double MAX_SINE_IMAG = 100.0; /** < Maximum value that is acceptable for cosh(z) doesn't overflow*/
 //Function declarations
void hypergeom(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag);

 //Function implementations
 /**
 * This small helper function computes the logarithm of the sine
 * Basically we ignore the small exponent of exp(i*(x+iy))=exp(-y+ix) and exp(-i*(x+iy))=exp(y-ix)
 * Depending on the sign of z_imag we can ignore one of these
 * */
void ln_sin_complex(double z_real,double z_imag,double* res_real,double* res_imag){
	if(z_imag>0){
		*res_real = z_imag-LOG_2;
		*res_imag = -z_real+MATH_PI_HALF;
		return;
	}else{
		*res_real = -z_imag-LOG_2;
		*res_imag = z_real-MATH_PI_HALF;
		return;
	}
}
 /**
 * This small helper function computes the logarithm of the cosine
 * Basically we ignore the small exponent of exp(i*(x+iy))=exp(-y+ix) and exp(-i*(x+iy))=exp(y-ix)
 * Depending on the sign of z_imag we can ignore one of these
 * */
void ln_cos_complex(double z_real,double z_imag, double* res_real, double* res_imag){
  if(z_imag>0){
    *res_real = z_imag-LOG_2;
    *res_imag = -z_real;
    return;
  }
  else{
    *res_real = -z_imag-LOG_2;
    *res_imag = z_real;
    return;
  }
}
/**
 * This function computes the real Gamma function, which we can simply take from the system library
 * 
 * */
double gamma_real(double z){
	return tgamma(z);
}
/**
 * This function computes the logarithm of the real Gamma function, which can be found in the system library
 * 
 * */
double ln_gamma_real(double z){
	return lgamma(z);
}
/**
 * This function computes the complex Gamma function using the Lanczos Approximation to order of NUM_GAMMA_COEFF, using the GAMMA_COEFF
 * The analytical continuation to Re(z)<0.5 is done using the reflection formula of the gamma function
 *
 * */
void gamma_complex(double z_real,double z_imag,double* res_real,double* res_imag){
	if (z_real > 0.5){
		z_real -= 1;
		//Initialize Lanczos sum
		double sum_real = 0.99999999999980993;
		double sum_imag = 0.0;
		double temp_abs_squared; double temp_real;
    int i;
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * temp_real / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t and log(t)
		double t_real = z_real + NUM_GAMMA_COEFF - 0.5;
		double t_abs_squared = t_real*t_real + z_imag*z_imag;
		double log_t_real = 0.5*log(t_abs_squared);
		double log_t_imag = atan2(z_imag, t_real);
		
		//Use t to calculate prefactor
		double exp_factor = exp((0.5 + z_real)*log_t_real - z_imag*log_t_imag - t_real);
		double temp = z_imag*log_t_real + (0.5 + z_real)*log_t_imag - z_imag;
		double exp_real = exp_factor*cos(temp);
		double exp_imag = exp_factor*sin(temp);
		
		//Calculate total
		*res_real = SQRT_2_PI*(exp_real*sum_real - exp_imag*sum_imag);
		*res_imag = SQRT_2_PI*(exp_imag*sum_real + exp_real*sum_imag);
		return;
	}
	else{
		//Reflection formula setup
		double sin_real = sin(MATH_PI*z_real)*cosh(MATH_PI*z_imag);
		double sin_imag = cos(MATH_PI*z_real)*sinh(MATH_PI*z_imag);
		z_real = -z_real;
		z_imag = -z_imag;
		
		//Initialize Lanczos sum
		double sum_real = 0.99999999999980993;
		double sum_imag = 0.0;
		double temp_abs_squared; double temp_real;
    int i;
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * (z_real + i + 1) / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t and log(t)
		double t_real = z_real + NUM_GAMMA_COEFF - 0.5;
		double t_abs_squared = t_real*t_real + z_imag*z_imag;
		double log_t_real = 0.5*log(t_abs_squared);
		double log_t_imag = atan2(z_imag, t_real);
		
		//Use t to calculate prefactor
		double exp_factor = exp((0.5 + z_real)*log_t_real - z_imag*log_t_imag - t_real);
		double temp = z_imag*log_t_real + (0.5 + z_real)*log_t_imag - z_imag;
		double exp_real = exp_factor*cos(temp);
		double exp_imag = exp_factor*sin(temp);
		
		//Calculate preliminary result of Gamma(1-z)
		double y_real = SQRT_2_PI*(exp_real*sum_real - exp_imag*sum_imag);
		double y_imag = SQRT_2_PI*(exp_imag*sum_real + exp_real*sum_imag);
		
		//Reflect using sin(pi*z)
		double ysin_real = y_real*sin_real - y_imag*sin_imag;
		double ysin_imag = y_imag*sin_real + y_real*sin_imag;
		double ysin_abs_squared = ysin_real*ysin_real + ysin_imag*ysin_imag;
		
		//Return final result
		*res_real = MATH_PI / (ysin_abs_squared)*ysin_real;
		*res_imag = -MATH_PI / (ysin_abs_squared)*ysin_imag;
		return;
	}
}
/**
 * This function computes the logarithm of the complex gamma function using a simplified algorithm
 * The analytical continuation to Re(z)<0.5 is done using the reflection formula of the gamma function
 *
 *  0.5*ln(2pi) + (z-0.5)*ln(z)-z+z/2*ln(z*sinh(1/z)+1/810z^6)
 * where we removed the computationally expensive error checking term in 1/810z^6
 *
 * Turns out to be a bit slower than the Lanczos approximation
 * */
void ln_gamma_complex_sinh(double z_real, double z_imag, double* res_real, double* res_imag){
	if(z_real>0.5){
		//Calculate log(z)
		double z_abs_squared = z_real*z_real + z_imag*z_imag;
		double ln_z_real = 0.5*log(z_abs_squared);
		double ln_z_imag = atan2(z_imag, z_real);
		
		//Calculate first terms in sum
		double first_real = ln_z_real*(z_real - 0.5) - ln_z_imag*z_imag;
		double first_imag = ln_z_real*z_imag + ln_z_imag*(z_real - 0.5);
		//Calculate 1/z and then sinh(1/z)
		double div_z_real = z_real / z_abs_squared; 
		double div_z_imag = -z_imag / z_abs_squared;
		double sinh_real = sinh(div_z_real)*cos(div_z_imag); 
		double sinh_imag = cosh(div_z_real)*sin(div_z_imag);
		
		//Calculate z*sinh(1/z) and then log(z*sinh(1/z))
		double z_sinh_real = z_real*sinh_real - z_imag*sinh_imag;
		double z_sinh_imag = z_imag*sinh_real + z_real*sinh_imag;
		double z_sinh_abs_squared = z_sinh_real*z_sinh_real + z_sinh_imag*z_sinh_imag;
		double ln_z_sinh_real = 0.5*log(z_sinh_abs_squared);
		double ln_z_sinh_imag = atan2(z_sinh_imag, z_sinh_real);
		
		//Finally multiply with z/2 and give back result
		double second_real = (z_real*ln_z_sinh_real - z_imag*ln_z_sinh_imag)*0.5;
		double second_imag = (z_real*ln_z_sinh_imag + z_imag*ln_z_sinh_real)*0.5;
		*res_real = ln2PI_05 + first_real - z_real + second_real;
		*res_imag = first_imag - z_imag + second_imag;
		return;
	}else{
		//Calculate log(sin(z)) for reflection using either sin(pi*z) and then log, or ln_sin_complex function defined above
		double log_sin_real; double log_sin_imag;
		if(z_imag>MAX_SINE_IMAG || z_imag<-MAX_SINE_IMAG){
			ln_sin_complex(MATH_PI*z_real,MATH_PI*z_imag,&log_sin_real,&log_sin_imag);
		}else{
			double sin_real = sin(MATH_PI*z_real)*cosh(MATH_PI*z_imag);
			double sin_imag = cos(MATH_PI*z_real)*sinh(MATH_PI*z_imag);
			log_sin_real = 0.5*log(sin_real*sin_real+sin_imag*sin_imag);
			log_sin_imag = atan2(sin_imag,sin_real);
		}
		//Caclulate reflected z and log(z)
		z_real = 1.0-z_real;
		z_imag = -z_imag;
		double z_abs_squared = z_real*z_real + z_imag*z_imag;
		double ln_z_real = 0.5*log(z_abs_squared);
		double ln_z_imag = atan2(z_imag, z_real);
		
		//Calculate first terms
		double first_real = ln_z_real*(z_real - 0.5) - ln_z_imag*z_imag;
		double first_imag = ln_z_real*z_imag + ln_z_imag*(z_real - 0.5);
		
		//Calculate 1/z and then sinh(1/z)
		//Overflow is not expected since div_z_real should be a small quantity
		double div_z_real = z_real / z_abs_squared; 
		double div_z_imag = -z_imag / z_abs_squared;
		double sinh_real = sinh(div_z_real)*cos(div_z_imag); 
		double sinh_imag = cosh(div_z_real)*sin(div_z_imag);
		//Caclulate z*sinh(1/z) and then log(z*sinh(1/z))
		double z_sinh_real = z_real*sinh_real - z_imag*sinh_imag;
		double z_sinh_imag = z_imag*sinh_real + z_real*sinh_imag;
		double z_sinh_abs_squared = z_sinh_real*z_sinh_real + z_sinh_imag*z_sinh_imag;
		double ln_z_sinh_real = 0.5*log(z_sinh_abs_squared);
		double ln_z_sinh_imag = atan2(z_sinh_imag, z_sinh_real);
		//Multiply it by z/2
		double second_real = (z_real*ln_z_sinh_real - z_imag*ln_z_sinh_imag)*0.5;
		double second_imag = (z_real*ln_z_sinh_imag + z_imag*ln_z_sinh_real)*0.5;
		double pref_real = first_real - z_real + second_real;
		double pref_imag = first_imag - z_imag + second_imag;
		//Return the result
		*res_real = lnPI_HALF_05-log_sin_real-pref_real;
		*res_imag = -log_sin_imag-pref_imag;
		return;
	}
}

/**
 * This function computes the logarithm of the complex gamma function using the Lanczos Approximation to order of NUM_GAMMA_COEFF, using the GAMMA_COEFF
 * The analytical continuation to Re(z)<0.5 is done using the reflection formula of the gamma function
 *   Turns out to be a bit faster than the sinh(1/z) approximation
 * */
void ln_gamma_complex_2(double z_real,double z_imag,double* res_real,double* res_imag){
	if (z_real > 0.5){
		z_real -= 1;
		//Initialize Lanczos sum
		double sum_real = 0.99999999999980993;
		double sum_imag = 0.0;
		double temp_abs_squared; double temp_real;
    int i;
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * temp_real / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t, log(t) and log(sum)
		double t_real = z_real + NUM_GAMMA_COEFF - 0.5;
		double log_sum_real = 0.5*log(sum_real*sum_real+sum_imag*sum_imag);
		double log_sum_imag = atan2(sum_imag,sum_real);
		double log_t_real = 0.5*log(t_real*t_real+z_imag*z_imag);
		double log_t_imag = atan2(z_imag,t_real);
		//Return result
		*res_real = ln2PI_05+log_sum_real-t_real+log_t_real*(0.5+z_real)-log_t_imag*z_imag;
		*res_imag = log_sum_imag-z_imag+log_t_imag*(0.5+z_real)+log_t_real*z_imag;
		return;
	}
	else{
		//Reflection formula sin setup
		double log_sin_real; double log_sin_imag;
		if(z_imag>MAX_SINE_IMAG || z_imag<-MAX_SINE_IMAG){
			ln_sin_complex(MATH_PI*z_real,MATH_PI*z_imag,&log_sin_real,&log_sin_imag);
		}else{
			double sin_real = sin(MATH_PI*z_real)*cosh(MATH_PI*z_imag);
			double sin_imag = cos(MATH_PI*z_real)*sinh(MATH_PI*z_imag);
			log_sin_real = 0.5*log(sin_real*sin_real+sin_imag*sin_imag);
			log_sin_imag = atan2(sin_imag,sin_real);
		}	
		//Calculate reflected z
		z_real = -z_real;
		z_imag = -z_imag;
		//Initialize Lanczos formula
		double sum_real = 0.99999999999980993;
		double sum_imag = 0.0;
		double temp_abs_squared; double temp_real;
    int i;
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * temp_real / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t, log(t) and log(sum)
		double t_real = z_real + NUM_GAMMA_COEFF - 0.5;
		double log_sum_real = 0.5*log(sum_real*sum_real+sum_imag*sum_imag);
		double log_sum_imag = atan2(sum_imag,sum_real);
		double log_t_real = 0.5*log(t_real*t_real+z_imag*z_imag);
		double log_t_imag = atan2(z_imag,t_real);
		//Return result
		*res_real = lnPI_HALF_05-log_sum_real+t_real-log_t_real*(0.5+z_real)+log_t_imag*z_imag-log_sin_real;
		*res_imag = -log_sum_imag+z_imag-log_t_imag*(0.5+z_real)-log_t_real*z_imag-log_sin_imag;
		return;
	}
}

 /**
 * This function calculates the hpyergeometric gaussian function for a,b complex, c,z real
 * This function uses the simple series representation
 * */
void hypergeom_series_sum_ab(double a_real, double a_imag, double b_real, double b_imag, double c, double z, double* result_real, double* result_imag){
	double C_real = 1.0; double C_imag = 0.0; //< Current term in sum
	double S_real = 1.0; double S_imag = 0.0; //< Current sum accumulation
	int i = 0;
	double temp_real = C_real;
	double temp_imag = C_imag;
	double ab_real,ab_imag;
	double C_factor;
  //printf("a= %.10e+%.10ej , b = %.10e+%.10ej , c = %.10e, z = %.10e\n",a_real,a_imag,b_real,b_imag,c,z);
	//Stop when convergence or maximum iterations reached
	while ((C_real*C_real + C_imag*C_imag) > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY*(S_real*S_real + S_imag*S_imag)  && i <= HYPERGEOM_MAX_ITER){ 
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		C_factor = z / (i + 1) / (c + i);
		//Result
		C_real = temp_real*C_factor;
		C_imag = temp_imag*C_factor;
		S_real += C_real;
		S_imag += C_imag;
		i += 1;
    //printf("C = %.10e+%.10ej and S = %.10e+%.10ej \n",C_real,C_imag,S_real,S_imag);
    
	}
  //printf("ab %4d values used Sum = %.20e+%.20ej\n",i,S_real,S_imag);
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * This function uses the simple series representation
 * */
 #include <stdio.h>
void hypergeom_series_sum_abc(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag){
	double C_real = 1.0; double C_imag = 0.0; /**< Current term in sum*/
	double S_real = 1.0; double S_imag = 0.0; /**< Current sum accumulation*/
	int i = 0;
	double temp_real = C_real;
	double temp_imag = C_imag;
	double ab_real,ab_imag;
	double C_factor;
	double den_abs_squared;
	while ((C_real*C_real + C_imag*C_imag) > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY *(S_real*S_real + S_imag*S_imag) && i <= HYPERGEOM_MAX_ITER){ //Needs more stopping terms
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		den_abs_squared = (c_real + i)*(c_real + i) + c_imag*c_imag;
		C_factor = z / (i + 1) / den_abs_squared;
		//Result
		C_real = (temp_real*(c_real + i) + temp_imag*c_imag)*C_factor;
		C_imag = (temp_imag*(c_real + i) - temp_real*c_imag)*C_factor;
		S_real += C_real;
		S_imag += C_imag;
		i += 1;
	}
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}
 /**
 * This function calculates the hpyergeometric gaussian function for a,b complex, c,z real
 * This function uses the simple series representation
 * */
void hypergeom_series_sum_ab_safe(double a_real, double a_imag, double b_real, double b_imag, double c, double z, double* result_real, double* result_imag, int* overflows){
	double C_real = 1.0; double C_imag = 0.0; //< Current term in sum
	double S_real = 1.0; double S_imag = 0.0; //< Current sum accumulation
	int i = 0;
  int count = 0;
  *overflows = 0;
	double temp_real = C_real;
	double temp_imag = C_imag;
	double ab_real,ab_imag;
	double C_factor;
  double C_abs = 1.0,S_abs=1.0;
  //double HYPERGEOM_SKIP_CRITERION = 0.01;
  //double HYPERGEOM_SUM_STEPS_SKIP = 4;
  //printf("a= %.10e+%.10ej , b = %.10e+%.10ej , c = %.10e, z = %.10e\n",a_real,a_imag,b_real,b_imag,c,z);
	//Stop when convergence or maximum iterations reached
	while (C_abs > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY*S_abs && count <= HYPERGEOM_MAX_ITER){ 
		//Calculate numerator
    ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
    ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
    temp_real = C_real*ab_real - C_imag*ab_imag;
    temp_imag = C_imag*ab_real + C_real*ab_imag;
    //Calculate denominator
    C_factor = z / (i + 1) / (c + i);
    //Result
    //if(i>9000 && c==1001.5 && fabs(C_real-temp_real*C_factor)<fabs(C_real)*HYPERGEOM_SKIP_CRITERION){
      /*double delta_real = temp_real*C_factor/C_real;
      double delta_imag = 0.0;//temp_imag*C_factor/C_imag;
      C_real = temp_real*C_factor;
      C_imag = temp_imag*C_factor;
      //printf("delta = %.10e , Predicted : %.10e , gotten %.10e \n",delta_real,delta_real*C_real,C_real);
      double temp_real = C_real;
      double temp_imag = C_imag;
      for(int j=0;j<HYPERGEOM_SUM_STEPS_SKIP;++j){
        S_real += C_real;
        S_imag += C_imag;
        C_real *=delta_real;
        C_imag *=delta_imag;
        i += 1;
      }*/
      //C_real = temp_real;
      //C_imag = temp_imag;
    //}
    //else{
      C_real = temp_real*C_factor;
      C_imag = temp_imag*C_factor;
    //}
    C_abs = C_real*C_real+C_imag*C_imag;
    S_abs = S_real*S_real+S_imag*S_imag;
    if(C_abs>ABS_LIMIT*ABS_LIMIT || S_abs > ABS_LIMIT*ABS_LIMIT){
      C_real/=ABS_LIMIT;
      C_imag/=ABS_LIMIT;
      S_real/=ABS_LIMIT;
      S_imag/=ABS_LIMIT;
      (*overflows)+=1;
      //printf("Overflow detected, but everything is fine \n");
    }
    S_real += C_real;
    S_imag += C_imag;
    i+=1;
    count+=1;
    
	}
  //printf("ab %4d values used Sum = %.20e+%.20ej\n",i,S_real,S_imag);
  
  *result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * This function uses the simple series representation
 * */
void hypergeom_series_sum_abc_safe(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag,int* overflows){
	double C_real = 1.0; double C_imag = 0.0; /**< Current term in sum*/
	double S_real = 1.0; double S_imag = 0.0; /**< Current sum accumulation*/
	int i = 0;
  *overflows = 0;
	double temp_real = C_real;
	double temp_imag = C_imag;
	double ab_real,ab_imag;
	double C_factor;
  double C_abs = 1.0,S_abs=1.0;
	double den_abs_squared;
	while (C_abs > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY *S_abs && i <= HYPERGEOM_MAX_ITER){ //Needs more stopping terms
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		den_abs_squared = (c_real + i)*(c_real + i) + c_imag*c_imag;
		C_factor = z / (i + 1) / den_abs_squared;
		//Result
		C_real = (temp_real*(c_real + i) + temp_imag*c_imag)*C_factor;
		C_imag = (temp_imag*(c_real + i) - temp_real*c_imag)*C_factor;
		S_real += C_real;
		S_imag += C_imag;
    if(C_abs>ABS_LIMIT*ABS_LIMIT || S_abs > ABS_LIMIT*ABS_LIMIT){
      C_real/=ABS_LIMIT;
      C_imag/=ABS_LIMIT;
      S_real/=ABS_LIMIT;
      S_imag/=ABS_LIMIT;
      (*overflows)+=1;
    }
		i += 1;
	}
  //printf("abc %4d values used Sum = %.10e+%.10ej\n",i,S_real,S_imag);
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * This function uses the Gosper algorithm to calculate it (Hypergeo)
 * */
void hypergeom_series_gosper(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* res_real, double* res_imag){
	//The three variables d,e,f of the gosper algorithm
	double d_real = 0.0; double d_imag = 0.0; 
	double e_real = 1.0; double e_imag = 0.0; 
	double f_real = 0.0; double f_imag = 0.0;
	int i = 0;
	//The next values of each variable
	double dnew_real,dnew_imag;
	double enew_real,enew_imag;
	double fnew_real,fnew_imag;
	
	while (i <= HYPERGEOM_MAX_ITER){
		//Recurrent factor
		double z_fact = z/(1-z);
		
		//Calculate the two denominators and their absolute values
		double den1_real = (i+1)*((2*i+c_real)*(2*i+c_real+1)-c_imag*c_imag);
		double den1_imag = (i+1)*((2*i+c_real)*c_imag+c_imag*(2*i+c_real+1));
		double den2_real = (2*i+c_real)*(1-z);
		double den2_imag = c_imag*(1-z);
		double den1_abs_squared = den1_real*den1_real+den1_imag*den1_imag;
		double den2_abs_squared = den2_real*den2_real+den2_imag*den2_imag;
		
		//Calculate dnew
		double dfact_real = e_real-((i+c_real-b_real-a_real)*d_real-(c_imag-b_imag-a_imag)*d_imag)*z_fact;
		double dfact_imag = e_imag-((i+c_real-b_real-a_real)*d_imag+(c_imag-b_imag-a_imag)*d_real)*z_fact;
		double abfact_real = (i+a_real)*(i+b_real)-a_imag*b_imag;
		double abfact_imag = (i+a_real)*b_imag+a_imag*(i+b_real);
		double dnum_real = (abfact_real*dfact_real-abfact_imag*dfact_imag)*z;
		double dnum_imag = (abfact_real*dfact_imag+abfact_imag*dfact_real)*z;
		dnew_real = (dnum_real*den1_real+dnum_imag*den1_imag)/den1_abs_squared;
		dnew_imag = (dnum_imag*den1_real-dnum_real*den1_imag)/den1_abs_squared;
		
		//Calculate enew
		double abd_fact_real = (a_real*b_real*d_real-a_real*b_imag*d_imag-a_imag*b_imag*d_real-a_imag*b_real*d_imag)*z_fact;
		double abd_fact_imag = (a_real*b_imag*d_real+a_imag*b_real*d_real+a_real*b_real*d_imag-a_imag*b_imag*d_imag)*z_fact;
		double efact_real = abd_fact_real+(i+c_real)*e_real-c_imag*e_imag;
		double efact_imag = abd_fact_imag+(i+c_real)*e_imag+c_imag*e_real;
		double enum_real = (abfact_real*efact_real-abfact_imag*efact_imag)*z;
		double enum_imag = (abfact_real*efact_imag+abfact_imag*efact_real)*z;
		enew_real = (enum_real*den1_real+enum_imag*den1_imag)/den1_abs_squared;
		enew_imag = (enum_imag*den1_real-enum_real*den1_imag)/den1_abs_squared;
		
		//Calculate fnew
		double first_real = i*((c_real-b_real-a_real)*z+i*(z-2)-c_real)-(a_real*b_real-a_imag*b_imag)*z;
		double first_imag = i*((c_imag-b_imag-a_imag)*z-c_imag)-(a_imag*b_real+a_real*b_imag)*z;
		double fnum_real = d_real*first_real-d_imag*first_imag;
		double fnum_imag = d_real*first_imag+d_imag*first_real;
		double fadd_real = (fnum_real*den2_real+fnum_imag*den2_imag)/den2_abs_squared;
		double fadd_imag = (fnum_imag*den2_real-fnum_real*den2_imag)/den2_abs_squared;
		fnew_real = f_real-fadd_real+e_real;
		fnew_imag = f_imag-fadd_imag+e_imag;
		
		//Compare fnew to f for convergence
		double fnew_abs = fnew_real*fnew_real+fnew_imag*fnew_imag;
		if(((fnew_real-f_real)*(fnew_real-f_real)+(fnew_imag-f_imag)*(fnew_imag-f_imag))<TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY*fnew_abs){
			*res_real = fnew_real; *res_imag = fnew_imag;
			return;
		}
		//Assign new variables
		d_real = dnew_real;
		d_imag = dnew_imag;
		e_real = enew_real;
		e_imag = enew_imag;
		f_real = fnew_real;
		f_imag = fnew_imag;
		i += 1;
	}
	*res_real = f_real;
	*res_imag = f_imag;
	return;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * This function uses a fractional representation
 * */
void hypergeom_series_fractional(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* res_real, double* res_imag){
	double alphabeta_real=1.0,alphabeta_imag=0.0;
	double beta_real=1.0,beta_imag=0.0;
	double gamma_real=1.0,gamma_imag=0.0;
	double temp;
	double cj_real,cj_imag;
	double ab_real,ab_imag;
	double curres_real,curres_imag;
	double preres_real=1.0,preres_imag=0.0;
	double gamma_abs_squared=1.0;
	double deltacur_real,deltacur_imag;
  int j;
	for(j=0;j<HYPERGEOM_MAX_ITER;++j){
		//Calculate the current denominator
		cj_real = (j+1)*(c_real+j);
		cj_imag = (j+1)*c_imag;
		//Calculate the current numerator
		ab_real = ((a_real+j)*(b_real+j)-a_imag*b_imag)*z;
		ab_imag = ((a_real+j)*b_imag+a_imag*(b_real+j))*z;
		//Calculate the total denominator and its absolute value
		temp = gamma_real;
		gamma_real = temp*cj_real-gamma_imag*cj_imag;
		gamma_imag = temp*cj_imag+gamma_imag*cj_real;
		gamma_abs_squared = gamma_real*gamma_real+gamma_imag*gamma_imag;
		//Calculate the new partial numerator beta
		temp = beta_real;
		beta_real = ab_real*temp-ab_imag*beta_imag;
		beta_imag = ab_real*beta_imag+ab_imag*temp;
		
		//Calculate the new numerator (uses beta calculated above, do not move)
		temp = alphabeta_real;
		alphabeta_real = temp*cj_real-alphabeta_imag*cj_imag+beta_real;
		alphabeta_imag = temp*cj_imag+alphabeta_imag*cj_real+beta_imag;
		
		//Check for necessity of resizing due to double overflow
		if(alphabeta_real*alphabeta_real+alphabeta_imag*alphabeta_imag>ABS_LIMIT || gamma_abs_squared>ABS_LIMIT){
			alphabeta_real/=ABS_LIMIT;
			alphabeta_imag/=ABS_LIMIT;
			beta_real/=ABS_LIMIT;
			beta_imag/=ABS_LIMIT;
			gamma_real/=ABS_LIMIT;
			gamma_imag/=ABS_LIMIT;
			gamma_abs_squared/=ABS_LIMIT*ABS_LIMIT;
		}
		
		//Calculate the current result
		curres_real = (alphabeta_real*gamma_real+alphabeta_imag*gamma_imag)/gamma_abs_squared;
		curres_imag = (alphabeta_imag*gamma_real-alphabeta_real*gamma_imag)/gamma_abs_squared;
		//Check for convergence of current result
		deltacur_real = curres_real-preres_real;
		deltacur_imag = curres_imag-preres_imag;
		if(deltacur_real*deltacur_real+deltacur_imag*deltacur_imag<HYPERGEOM_EPSILON*(preres_real*preres_real+preres_imag*preres_imag)){
			break;
		}
		//Assign new variables (do not move)
		preres_real = curres_real;
		preres_imag = curres_imag;
	}
	*res_real = curres_real;
	*res_imag = curres_imag;
	return;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * This function uses an expansion around z = 0.5
 * */
void hypergeom_series_zexpansion(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* res_real, double* res_imag){
	//Recurrent factors
	double z_fac = (z/(z-2.0));
	double c_abs_squared = c_real*c_real+c_imag*c_imag;
	//Declare temporary variables
	double c_cur_real,c_cur_abs_squared;
	double phi_prev_real = 1.0; 
	double phi_prev_imag=0.0;
	
	double phi_cur_real = 1.0-2.0*(b_real*c_real+b_imag*c_imag)/c_abs_squared;
	double phi_cur_imag = -2.0*(b_imag*c_real-b_real*c_imag)/c_abs_squared;
	
	double aterm_real = a_real*z_fac;
	double aterm_imag = a_imag*z_fac;
	
	double phi_next_real,phi_next_imag;
	double temp_real,temp_imag;
	
	double sum_real = 1.0+(aterm_real*phi_cur_real-aterm_imag*phi_cur_imag);
	double sum_imag =     (aterm_real*phi_cur_imag+aterm_imag*phi_cur_real);
	double cur_real = 1.0;
	double cur_imag = 0.0;
	int i;
	for(i=1;i<HYPERGEOM_MAX_ITER;++i){
		//Calculate new Gamma(a+n)/Gamma(a)*(z/(z-2))^n/n! term
		temp_real = aterm_real;
		aterm_real = (temp_real*(a_real+i)-aterm_imag*a_imag)/(i+1)*z_fac;
		aterm_imag = (temp_real*a_imag+(a_real+i)*aterm_imag)/(i+1)*z_fac;
		
		//Calculate new phi by using recursion relations
		temp_real = (i*phi_prev_real-(2*b_real-c_real)*phi_cur_real+(2*b_imag-c_imag)*phi_cur_imag);
		temp_imag = (i*phi_prev_imag-(2*b_real-c_real)*phi_cur_imag-(2*b_imag-c_imag)*phi_cur_real);
		c_cur_real = c_real+i;
		c_cur_abs_squared = c_cur_real*c_cur_real+c_imag*c_imag;
		phi_next_real = (temp_real*c_cur_real+temp_imag*c_imag)/c_cur_abs_squared;
		phi_next_imag = (temp_imag*c_cur_real-temp_real*c_imag)/c_cur_abs_squared;

		//Multiply both terms to obtain next term in sum
		cur_real = (aterm_real*phi_next_real-aterm_imag*phi_next_imag);
		cur_imag = (aterm_imag*phi_next_real+aterm_real*phi_next_imag);
		
		sum_real += cur_real;
		sum_imag += cur_imag;
		if((cur_real*cur_real+cur_imag*cur_imag)<HYPERGEOM_EPSILON*(sum_real*sum_real+sum_imag*sum_imag)){
			break;
		}
		//Assign new variables
		phi_prev_real = phi_cur_real;
		phi_prev_imag = phi_cur_imag;
		phi_cur_real = phi_next_real;
		phi_cur_imag = phi_next_imag;
		//Check for necessity of resizing
		if(phi_cur_real*phi_cur_real+phi_cur_imag*phi_cur_imag>ABS_LIMIT || aterm_real*aterm_real+aterm_imag*aterm_imag<INV_ABS_LIMIT){
			phi_cur_real*=INV_ABS_LIMIT;
			phi_cur_imag*=INV_ABS_LIMIT;
			phi_prev_real*=INV_ABS_LIMIT;
			phi_prev_imag*=INV_ABS_LIMIT;
			aterm_real*=ABS_LIMIT;
			aterm_imag*=ABS_LIMIT;
		}
	}
	//Calculate final prefactor of (1-0.5z)^(-a)
	double z_fac2 = log(1-0.5*z);
	double exp_factor = exp(-z_fac2*a_real);
	double exp_real = exp_factor*cos(z_fac2*a_imag);
	double exp_imag = -exp_factor*sin(z_fac2*a_imag);
	
	*res_real = sum_real*exp_real-sum_imag*exp_imag;
	*res_imag = sum_imag*exp_real+sum_real*exp_imag;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * This function uses recursion relations, making it do much more work
 * */
void hypergeom_recursionbc(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag){
	//Calculate the number of steps needed to get back to reasonable ranges B_C_MAX;
	int n_steps = ((int)(b_real)) - B_C_MAX + 1; 
	if(n_steps<=0){
			// TODO : ERROR MESSAGE
			// return __FAILURE__;
		}
	double b_start = b_real - n_steps - 1;
	double c_start = c_real - n_steps - 1;
	
	double cur_hyper_real = 0.0; double cur_hyper_imag = 0.0;
	hypergeom(a_real, a_imag, b_start+1, b_imag, c_start+1, c_imag, z, &cur_hyper_real, &cur_hyper_imag);
	double prev_hyper_real = 0.0; double prev_hyper_imag = 0.0;
	hypergeom(a_real, a_imag, b_start, b_imag, c_start, c_imag, z, &prev_hyper_real, &prev_hyper_imag);
	double F_previous_real = prev_hyper_real; double F_previous_imag = prev_hyper_imag;
	double F_current_real = cur_hyper_real; double F_current_imag = cur_hyper_imag;
	int i;
  for (i = 0; i < n_steps; ++i){
		double c_cur = c_start + i;
		double b_cur = b_start + i;
		
		double pref_num_real = (c_cur)*(c_cur + 1) - c_imag*c_imag;
		double pref_num_imag = 2 * (c_cur)*c_imag + c_imag;
		
		double pref_den_real = ((b_cur + 1)*(c_cur - a_real + 1) - b_imag*(c_imag - a_imag));
		double pref_den_imag = (b_imag*(c_cur - a_real + 1) + (b_cur + 1)*(c_imag - a_imag));
		
		double pref_den_abs_squared = pref_den_real*pref_den_real + pref_den_imag*pref_den_imag;
		double pref_real = (pref_num_real*pref_den_real + pref_num_imag*pref_den_imag) / pref_den_abs_squared / z;
		double pref_imag = (pref_num_imag*pref_den_real - pref_num_real*pref_den_imag) / pref_den_abs_squared / z;
		
		double c_abs_squared = c_cur*c_cur + c_imag*c_imag;
		
		double pref_1_real = 1.0 - ((a_real - b_cur - 1)*c_cur + (a_imag - b_imag)*c_imag) * z / c_abs_squared;
		double pref_1_imag = ((a_real - b_cur - 1)*c_imag - (a_imag - b_imag)*c_cur)*z / c_abs_squared;
		
		double bracket_real = pref_1_real*F_current_real - pref_1_imag*F_current_imag - F_previous_real;
		double bracket_imag = pref_1_imag*F_current_real + pref_1_real*F_current_imag - F_previous_imag;
		
		double F_next_real = bracket_real*pref_real - bracket_imag*pref_imag;
		double F_next_imag = bracket_real*pref_imag + bracket_imag*pref_real;
		
		F_previous_real = F_current_real;
		F_previous_imag = F_current_imag;
		F_current_real = F_next_real;
		F_current_imag = F_next_imag;
	}
	*result_real = F_current_real;
	*result_imag = F_current_imag;
	return;
}
void hypergeom_recursionab(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag){
	int n_steps = ((int)(a_real))-A_B_MAX + 1;
	if(n_steps<=0){
		// TODO : ERROR MESSAGE
		// return __FAILURE__;
	}
	double a_start = a_real - n_steps - 1;
	double b_start = b_real + n_steps + 1;
	double cur_hyper_real = 0.0; double cur_hyper_imag = 0.0;
	hypergeom(a_start+1, a_imag, b_start-1, b_imag, c_real, c_imag, z, &cur_hyper_real, &cur_hyper_imag);
	double prev_hyper_real = 0.0; double prev_hyper_imag = 0.0;
	hypergeom(a_start, a_imag, b_start, b_imag, c_real, c_imag, z, &prev_hyper_real, &prev_hyper_imag);
	double F_previous_real = prev_hyper_real; double F_previous_imag = prev_hyper_imag;
	double F_current_real = cur_hyper_real; double F_current_imag = cur_hyper_imag;
	int i;
  for (i = 0; i < n_steps; ++i){
		double a_cur = a_start + i + 1.0;
		double b_cur = b_start - i - 1.0;
		
		double ba1_real = b_cur - a_cur - 1; double ba1_imag = b_imag - a_imag;
		double ab1_real = b_cur - a_cur + 1; double ab1_imag = b_imag - a_imag;
		double ab1_abs_squared = (ab1_real*ab1_real + ab1_imag*ab1_imag);
		
		double gamma_real = (ba1_real*b_cur*ab1_real-ba1_imag*b_imag*ab1_real+ba1_imag*b_cur*ab1_imag+ba1_real*b_imag*ab1_imag)/ab1_abs_squared;
		double gamma_imag = (ba1_imag*b_cur*ab1_real+ba1_real*b_imag*ab1_real-ba1_real*b_cur*ab1_imag+ba1_imag*b_imag*ab1_imag)/ab1_abs_squared;

		double first_real = a_cur*(c_real-a_cur-1)-a_imag*(c_imag-a_imag);
		double first_imag = a_imag*(c_real-a_cur-1)+a_cur*(c_imag-a_imag);
		double second_real = (ba1_real*(a_cur-b_cur)-ba1_imag*(a_imag-b_imag))*(z-1);
		double second_imag = (ba1_imag*(a_cur-b_cur)+ba1_real*(a_imag-b_imag))*(z-1);
		double third_real = gamma_real*(c_real-b_cur-1)-gamma_imag*(c_imag-b_imag);
		double third_imag = gamma_real*(c_imag-b_imag)+gamma_imag*(c_real-b_cur-1);

		
		double alpha_real = first_real+second_real+third_real;
		double alpha_imag = first_imag+second_imag+third_imag;
		double beta_real = gamma_real*(a_cur-c_real)-gamma_imag*(a_imag-c_imag);
		double beta_imag = gamma_imag*(a_cur-c_real)+gamma_real*(a_imag-c_imag);
		double div_pref_real = a_cur*(c_real - b_cur)-a_imag*(c_imag-b_imag);
		double div_pref_imag = a_imag*(c_real - b_cur)+a_cur*(c_imag-b_imag);
		double div_pref_abs_squared = div_pref_real*div_pref_real + div_pref_imag*div_pref_imag;

		double tot_ap_bm_real = alpha_real*F_current_real - alpha_imag*F_current_imag;
		double tot_ab_real = beta_real*F_previous_real - beta_imag*F_previous_imag;
		double tot_ap_bm_imag = alpha_imag*F_current_real + alpha_real*F_current_imag;
		double tot_ab_imag = beta_imag*F_previous_real + beta_real*F_previous_imag;
		double F_next_real = ((tot_ap_bm_real + tot_ab_real)*div_pref_real + (tot_ap_bm_imag + tot_ab_imag)*div_pref_imag) / div_pref_abs_squared;
		double F_next_imag = ((tot_ap_bm_imag + tot_ab_imag)*div_pref_real - (tot_ap_bm_real + tot_ab_real)*div_pref_imag) / div_pref_abs_squared;
		F_previous_real = F_current_real;
		F_previous_imag = F_current_imag;
		F_current_real = F_next_real;
		F_current_imag = F_next_imag;
	}
	*result_real = F_current_real;
	*result_imag = F_current_imag;
	return;

}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * It selects from previously defined functions
 * 
 * For Re(a,b,c)>A_B_C_MAX, we use recursion relations, otherwise taylor expansions
 * 
 * There are only two cases of recursion relations we actually need to deal with in the course of this program
 * One from b,c ~ l , so we have b and c with large real part, thus using b,c recursion relations.
 * One from a,-b ~ l, so we have a and -b with large real part, thus using a,-b recursion relations
 *
 * As it currently stands, the b,c recursion relation is numerically unstable, since it depends on the exact difference between two large hypergeometric functions
 * 
 * This function is not guaranteed to give a reasonable result, especially for z->1. 
 * The handling of the convergence is instead done by the bessel_integral functions!
 * */

void hypergeom(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag){

	//Should the border value of z at 0 be used, we know the value of the hypergeometric function immediately
	if(z==0.0){
		*result_real = 1; *result_imag=0; return;
	}
	
	//b-c small, but b,c large (use recursion on b,c) 
	if (b_real > B_C_MAX && c_real > B_C_MAX){
		hypergeom_recursionbc(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
	}
	
	//a+b small, but a,-b large (use recursion on a,-b)
	if (a_real > A_B_MAX && b_real < A_B_MAX){
		hypergeom_recursionab(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
	}
	
	//If no reduction case is required, we use the simple taylor approximation
	//In the case of c_imag ==0.0, we should use the simplified version to save calculation time
	if(c_imag==0.0){
		hypergeom_series_sum_ab(a_real,a_imag,b_real,b_imag,c_real,z,result_real,result_imag);
		//hypergeom_series_gosper(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
    
    //hypergeom_series_sum_abc(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
		return;	
	}
	//In the general case, we use the general version of the taylor approximation
	else{
		hypergeom_series_sum_abc(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
		return;
	}
}
void hypergeom_safe(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag,int* overflows){
  //Should the border value of z at 0 be used, we know the value of the hypergeometric function immediately
	if(z==0.0){
		*result_real = 1; *result_imag=0; *overflows=0; return;
	}
	
	//If no reduction case is required, we use the simple taylor approximation
	//In the case of c_imag ==0.0, we should use the simplified version to save calculation time
	if(c_imag==0.0){
		hypergeom_series_sum_ab_safe(a_real,a_imag,b_real,b_imag,c_real,z,result_real,result_imag,overflows);
		//hypergeom_series_gosper(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
    
    //hypergeom_series_sum_abc(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
		return;	
	}
	//In the general case, we use the general version of the taylor approximation
	else{
		hypergeom_series_sum_abc_safe(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag,overflows);
		return;
	}
}

void hypergeom3F2_specific(double b_real,double c_real,double c_imag, double d_real,double d_imag, double e_real, double e_imag,double* res_real, double* res_imag){
  double cur_real = 1.0; double cur_imag = 0.0; /**< Current term in sum*/
	double S_real = 1.0; double S_imag = 0.0; /**< Current sum accumulation*/
	int i = 0;
  double ccur_real,ccur_imag;
  double de_real,de_imag;
	double cur_factor;
	double den_abs_squared;
  double num_real,num_imag;
	while (
      (cur_real*cur_real + cur_imag*cur_imag) > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY *(S_real*S_real + S_imag*S_imag)
       && i <= HYPERGEOM_MAX_ITER
    ){ 
		//Calculate numerator
    ccur_real = (cur_real*(c_real+i)-cur_imag*c_imag);
    ccur_imag = (cur_real*c_imag+cur_imag*(c_real+i));
    num_real = (b_real+i)*ccur_real;
    num_imag = (b_real+i)*ccur_imag;
    //Calculate denominator
    de_real = (d_real+i)*(e_real+i) - d_imag*e_imag;
    de_imag = (d_real+i)*e_imag + d_imag*(e_real+i);
		den_abs_squared = de_real*de_real+de_imag*de_imag;
		//Calculate prefactor
		cur_factor = 1.0 / den_abs_squared;
		//Result
		cur_real = (num_real*de_real + num_imag*de_imag)*cur_factor;
		cur_imag = (num_imag*de_real - num_real*de_imag)*cur_factor;
		S_real += cur_real;
		S_imag += cur_imag;
		i += 1;
  }
	*res_real = S_real;
	*res_imag = S_imag;
	return;
}
void bessel_integral_integral(double p, double l, double nu_real, double nu_imag,double* res_real,double* res_imag){
  double prefactor = MATH_PI_3_2;
  double pref_real,pref_imag;
  if(nu_imag<NU_IMAG_LIMIT){
    if(l<BESSEL_INTEGRAL_MAX_L){
      double gamma_num_first_real,gamma_num_first_imag;
      double gamma_num_scnd_real,gamma_num_scnd_imag;
      double gamma_den_first_real,gamma_den_first_imag;
      double gamma_den_scnd_real,gamma_den_scnd_imag;
      gamma_complex(2.0-nu_real*0.5,-nu_imag*0.5,&gamma_num_first_real,&gamma_num_first_imag);
      gamma_complex(l+nu_real*0.5,nu_imag*0.5,&gamma_num_scnd_real,&gamma_num_scnd_imag);
      gamma_complex(2.5-nu_real*0.5,-nu_imag*0.5,&gamma_den_first_real,&gamma_den_first_imag);
      gamma_complex(3.0+l-nu_real*0.5,-nu_imag*0.5,&gamma_den_scnd_real,&gamma_den_scnd_imag);
      double num_real = gamma_num_first_real*gamma_num_scnd_real-gamma_num_first_imag*gamma_num_scnd_imag;
      double num_imag = gamma_num_first_imag*gamma_num_scnd_real+gamma_num_first_real*gamma_num_scnd_imag;
      double den_real = gamma_den_first_real*gamma_den_scnd_real-gamma_den_first_imag*gamma_den_scnd_imag;
      double den_imag = gamma_den_first_imag*gamma_den_scnd_real+gamma_den_first_real*gamma_den_scnd_imag;
      double den_abs_squared = den_real*den_real+den_imag*den_imag;
      pref_real = (num_real*den_real+num_imag*den_imag)/den_abs_squared;
      pref_imag = (num_imag*den_real-num_real*den_imag)/den_abs_squared;
    }
    else{
      double gamma_num_first_real,gamma_num_first_imag;
      double gamma_den_first_real,gamma_den_first_imag;
      gamma_complex(2.0-nu_real*0.5,-nu_imag*0.5,&gamma_num_first_real,&gamma_num_first_imag);
      gamma_complex(2.5-nu_real*0.5,-nu_imag*0.5,&gamma_den_first_real,&gamma_den_first_imag);
      double a1,b1,a3,b3;
      ln_gamma_complex_2(l+0.5*nu_real,0.5*nu_imag,&a1,&b1);
      ln_gamma_complex_2(3.0+l-0.5*nu_real,-0.5*nu_imag,&a3,&b3);
      double exp_factor = exp(a1-a3);
      double fac1_real = cos(b1-b3)*exp_factor;
      double fac1_imag = sin(b1-b3)*exp_factor;
      double den_abs_squared = gamma_den_first_real*gamma_den_first_real+gamma_den_first_imag*gamma_den_first_imag;
      double fac2_real = (gamma_num_first_real*gamma_den_first_real+gamma_num_first_imag*gamma_den_first_imag)/den_abs_squared;
      double fac2_imag = (gamma_num_first_imag*gamma_den_first_real-gamma_num_first_real*gamma_den_first_imag)/den_abs_squared;
      pref_real = fac1_real*fac2_real-fac1_imag*fac2_imag;
      pref_imag = fac1_real*fac2_imag+fac1_imag*fac2_real;
    }
  }
  else{
    double a0,b0,a1,b1,a2,b2,a3,b3;
    ln_gamma_complex_2(2.0-nu_real*0.5,-nu_imag*0.5,&a0,&b0);
    ln_gamma_complex_2(l+0.5*nu_real,0.5*nu_imag,&a1,&b1);
    ln_gamma_complex_2(2.5-nu_real*0.5,-nu_imag*0.5,&a2,&b2);
    ln_gamma_complex_2(3.0+l-0.5*nu_real,-0.5*nu_imag,&a3,&b3);
    double exp_factor = exp(a0+a1 - a2-a3);
    pref_real = exp_factor*cos(b0+b1 - b2-b3);
    pref_imag = exp_factor*sin(b0+b1 - b2-b3);
  }
  double hyp_real,hyp_imag;
  hypergeom3F2_specific(1.0+(l-p)*0.5,3.0-nu_real,-nu_imag,3.0+l-nu_real*0.5,-nu_imag*0.5,2.5-nu_real*0.5,-nu_imag*0.5,&hyp_real,&hyp_imag);
  *res_real = prefactor*(pref_real*hyp_real-pref_imag*hyp_imag);
  *res_imag = prefactor*(pref_real*hyp_imag+pref_imag*hyp_real);
  return;
}
/**
 * This function calculates the analytical limit for the bessel_integral for t->1 (z=t²->1 as well)
 * 
 * */
void bessel_analytic_limit_atz1(double l, double nu_real, double nu_imag, double* res_real, double* res_imag){
	double exp_real; double exp_imag;
	//Gamma-overflow safe evaluation of Gamma(l+nu/2)/Gamma(l+2-nu/2)
	if(l<=BESSEL_INTEGRAL_MAX_L && nu_imag < NU_IMAG_LIMIT){
		double num_gamma_real = 0.0; double num_gamma_imag = 0.0;
		gamma_complex(l + nu_real*0.5, nu_imag*0.5,&num_gamma_real,&num_gamma_imag);
		double den_gamma_real = 0.0; double den_gamma_imag = 0.0;
		gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5,&den_gamma_real,&den_gamma_imag);
		double den_abs_squared = den_gamma_real*den_gamma_real+den_gamma_imag*den_gamma_imag;
		exp_real = (num_gamma_real*den_gamma_real+num_gamma_imag*den_gamma_imag)/den_abs_squared;
		exp_imag = (num_gamma_imag*den_gamma_real-num_gamma_real*den_gamma_imag)/den_abs_squared;
	}
	else{
		double first_ln_gamma_real = 0.0; double first_ln_gamma_imag = 0.0;
		ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &first_ln_gamma_real, &first_ln_gamma_imag);
		double scnd_ln_gamma_real = 0.0; double scnd_ln_gamma_imag = 0.0;
		ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &scnd_ln_gamma_real, &scnd_ln_gamma_imag);
		double exp_factor = exp(first_ln_gamma_real - scnd_ln_gamma_real);
		exp_real = exp_factor*cos(first_ln_gamma_imag - scnd_ln_gamma_imag);
		exp_imag = exp_factor*sin(first_ln_gamma_imag - scnd_ln_gamma_imag);
	}
	//Gamma-overflow safe evaluation of Gamma(1-nu/2)/Gamma((3-nu)/2)
	double div_real,div_imag;
	if(nu_imag < NU_IMAG_LIMIT){
		double first_gamma_real = 0.0; double first_gamma_imag = 0.0;
		gamma_complex(1-nu_real*0.5, -nu_imag*0.5, &first_gamma_real, &first_gamma_imag);
		double scnd_gamma_real = 0.0; double scnd_gamma_imag = 0.0;
		gamma_complex(1.5- nu_real*0.5, -nu_imag*0.5, &scnd_gamma_real, &scnd_gamma_imag);
		
		double scnd_gamma_abs_squared = scnd_gamma_real*scnd_gamma_real + scnd_gamma_imag*scnd_gamma_imag;
		div_real = (first_gamma_real*scnd_gamma_real + first_gamma_imag*scnd_gamma_imag) / scnd_gamma_abs_squared;
		div_imag = (first_gamma_imag*scnd_gamma_real - first_gamma_real*scnd_gamma_imag) / scnd_gamma_abs_squared;
	}
	else{
		double first_gamma_real = 0.0; double first_gamma_imag = 0.0;
		ln_gamma_complex_2(1-nu_real*0.5, -nu_imag*0.5, &first_gamma_real, &first_gamma_imag);
		double scnd_gamma_real = 0.0; double scnd_gamma_imag = 0.0;
		ln_gamma_complex_2(1.5- nu_real*0.5, -nu_imag*0.5, &scnd_gamma_real, &scnd_gamma_imag);
		
		double g_exp_factor = exp(first_gamma_real-scnd_gamma_real);
		double g_exp_phase = first_gamma_imag-scnd_gamma_imag;
		div_real = g_exp_factor*cos(g_exp_phase);
		div_imag = g_exp_factor*sin(g_exp_phase);
	}
	//Result
	*res_real = MATH_PI_3_2*(exp_real*div_real - exp_imag*div_imag);
	*res_imag = MATH_PI_3_2*(exp_imag*div_real + exp_real*div_imag);
	return;
}
/**
 * This function calculates the analytical limit for the bessel_integral of two different l for t->1 (z=t²->1 as well)
 * 
 * */
void double_bessel_analytic_limit_atz1(double l1, double l2, double nu_real, double nu_imag, double* res_real, double* res_imag){
	double exp_real; double exp_imag;
	//Gamma-overflow safe evaluation of Gamma(l+nu/2)/Gamma(l+2-nu/2)
	if(0.5*(l1+l2)<=BESSEL_INTEGRAL_MAX_L && nu_imag < NU_IMAG_LIMIT){
		double num_gamma_real = 0.0; double num_gamma_imag = 0.0;
		gamma_complex( (l1 + l2 + nu_real)*0.5, nu_imag*0.5,&num_gamma_real,&num_gamma_imag);
		double den_gamma_real = 0.0; double den_gamma_imag = 0.0;
		gamma_complex( (l1 + l2 + 4 - nu_real)*0.5, -nu_imag*0.5,&den_gamma_real,&den_gamma_imag);
		double den_abs_squared = den_gamma_real*den_gamma_real+den_gamma_imag*den_gamma_imag;
		exp_real = (num_gamma_real*den_gamma_real+num_gamma_imag*den_gamma_imag)/den_abs_squared;
		exp_imag = (num_gamma_imag*den_gamma_real-num_gamma_real*den_gamma_imag)/den_abs_squared;
	}
	else{
		double first_ln_gamma_real = 0.0; double first_ln_gamma_imag = 0.0;
		ln_gamma_complex_2( (l1 + l2 + nu_real)*0.5, nu_imag*0.5, &first_ln_gamma_real, &first_ln_gamma_imag);
		double scnd_ln_gamma_real = 0.0; double scnd_ln_gamma_imag = 0.0;
		ln_gamma_complex_2( (l1 + l2 + 4 - nu_real)*0.5, -nu_imag*0.5, &scnd_ln_gamma_real, &scnd_ln_gamma_imag);
		double exp_factor = exp(first_ln_gamma_real - scnd_ln_gamma_real);
		exp_real = exp_factor*cos(first_ln_gamma_imag - scnd_ln_gamma_imag);
		exp_imag = exp_factor*sin(first_ln_gamma_imag - scnd_ln_gamma_imag);
	}
	//Gamma-overflow safe evaluation of Gamma(1-nu/2)/Gamma((3-nu)/2)
	double div_real,div_imag;
	if(nu_imag < NU_IMAG_LIMIT){
		double first_gamma_real = 0.0; double first_gamma_imag = 0.0;
		gamma_complex(2-nu_real, -nu_imag, &first_gamma_real, &first_gamma_imag);
		double scnd_gamma_real = 0.0; double scnd_gamma_imag = 0.0;
		gamma_complex(1.5 - nu_real*0.5 + (l1-l2)*0.5, -nu_imag*0.5, &scnd_gamma_real, &scnd_gamma_imag);
    double third_gamma_real = 0.0; double third_gamma_imag = 0.0;
    gamma_complex(1.5 - nu_real*0.5 + (l2-l1)*0.5, -nu_imag*0.5, &third_gamma_real, &third_gamma_imag);
    double exp_factor = exp(LOG_2*(nu_real-1));
    double exp2_real = exp_factor*cos(LOG_2*nu_imag);
    double exp2_imag = exp_factor*sin(LOG_2*nu_imag);
    
    double num_real = (first_gamma_real*exp2_real-first_gamma_imag*exp2_imag);
    double num_imag = (first_gamma_real*exp2_imag+first_gamma_imag*exp2_real);
    
    double den_real = (scnd_gamma_real*third_gamma_real-scnd_gamma_imag*third_gamma_imag);
    double den_imag = (scnd_gamma_real*third_gamma_imag+first_gamma_imag*third_gamma_real);
    
    double den_abs_squared = den_real*den_real + den_imag*den_imag;
    div_real = (num_real*den_real+num_imag*den_imag)/den_abs_squared;
    div_imag = (num_imag*den_real-num_real*den_imag)/den_abs_squared;
	}
	else{
		double first_gamma_real = 0.0; double first_gamma_imag = 0.0;
		ln_gamma_complex_2(2-nu_real, -nu_imag, &first_gamma_real, &first_gamma_imag);
		double scnd_gamma_real = 0.0; double scnd_gamma_imag = 0.0;
		ln_gamma_complex_2(1.5 - nu_real*0.5 + (l1-l2)*0.5, -nu_imag*0.5, &scnd_gamma_real, &scnd_gamma_imag);
		double third_gamma_real = 0.0; double third_gamma_imag = 0.0;
    ln_gamma_complex_2(1.5 - nu_real*0.5 + (l2-l1)*0.5, -nu_imag*0.5, &third_gamma_real, &third_gamma_imag);
		
		double g_exp_factor = exp(first_gamma_real+LOG_2*(nu_real-1)-scnd_gamma_real-third_gamma_real);
		double g_exp_phase = first_gamma_imag+LOG_2*nu_imag-scnd_gamma_imag-third_gamma_imag;
		div_real = g_exp_factor*cos(g_exp_phase);
		div_imag = g_exp_factor*sin(g_exp_phase);
	}
	//Result
	*res_real = MATH_PI*MATH_PI*(exp_real*div_real - exp_imag*div_imag);
	*res_imag = MATH_PI*MATH_PI*(exp_imag*div_real + exp_real*div_imag);
	return;
}
/**
 * This function calculates the integral over the bessel functions given in equation B.1
 * We have to be careful about following facts: 
 *
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where (1-t*t)<<1 is not true
 * Special care has to be taken because of that, which is explained further in the paper. 
 * Here, the CALLING function will have care about always keeping (1-t*t)<<1 .
 *
 * This function is mostly applicable in the low t, low nu regime
 * (Using bessel_integral_lowt_transform should be preferred)
 * */
 void bessel_integral_lowt_lownu(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Check for large imaginary nu to safely calculate gamma functions
	if(nu_imag>NU_IMAG_LIMIT){
		
		//Calculate total prefactor of 2^(nu-1)
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex_2(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		
		//Divide result by overflow-safe Gamma(l+3/2)
		double exp_real,exp_imag;
		if(l < BESSEL_INTEGRAL_MAX_L){
			double scnd_den_gamma = gamma_real(l + 1.5);
			exp_factor = exp(lgamma_num_real - a_0+log(t)*l)/scnd_den_gamma;
			
			exp_real = exp_factor*cos(lgamma_num_imag - b_0);
			exp_imag = exp_factor*sin(lgamma_num_imag - b_0);
		}else{
			double scnd_den_gamma = ln_gamma_real(l + 1.5);
			exp_factor = exp(lgamma_num_real - a_0 - scnd_den_gamma+log(t)*l);
			
			exp_real = exp_factor*cos(lgamma_num_imag - b_0);
			exp_imag = exp_factor*sin(lgamma_num_imag - b_0);
		}
		//We now obtained the total prefactor of the hypergeometric function
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*exp_real - pref_imag*exp_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*exp_real + pref_real*exp_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.5 - 0.5, nu_imag*0.5, l + nu_real*0.5, nu_imag*0.5, l + 1.5 ,0.0, t*t, &hypergeom_real, &hypergeom_imag);
		
		//Result
		*res_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
	
	else if (l < BESSEL_INTEGRAL_MAX_L){
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate Gamma((3-nu)/2) and Gamma(l+3/2) as the denominator
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		double scnd_den_gamma = gamma_real(l + 1.5);
		
		//Calculate the numerator
		double gamma_num_real = 0.0; double gamma_num_imag = 0.0;
		gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &gamma_num_real, &gamma_num_imag);
		
		//Calculate fraction of numerator/denominator and multiply with prefactor
		double frac_real = (first_den_gamma_real*gamma_num_real+first_den_gamma_imag*gamma_num_imag)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double frac_imag = (first_den_gamma_real*gamma_num_imag-first_den_gamma_imag*gamma_num_real)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.5 - 0.5, nu_imag*0.5, l + nu_real*0.5, nu_imag*0.5, l + 1.5 ,0.0, t*t, &hypergeom_real, &hypergeom_imag);
		
		//Final prefactor of t^l, which is going to make the integral vanishingly small for high l and low t
		double t_l_pref = pow(t,l);
		*res_real = t_l_pref*(hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = t_l_pref*(hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;		
	}
	else{
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Only Gamma((3-nu)/2) can be calculated directly
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;

		//Calculate Gamma(l+nu/2)/Gamma(l+3/2)
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		exp_factor = exp(lgamma_num_real - scnd_den_gamma); 
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);
		
		//Calculate fraction of numerator and denominator
		double frac_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		double frac_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
		
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.5 - 0.5, nu_imag*0.5, l + nu_real*0.5, nu_imag*0.5, l + 1.5 ,0.0, t*t, &hypergeom_real, &hypergeom_imag);
		
		//Final prefactor of t^l, which is going to make the integral vanishingly small for high l and low t
		double t_l_pref = pow(t,l);
		*res_real = t_l_pref*(hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = t_l_pref*(hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
}
/**
 * This function calculates the integral over the bessel functions given in equation B.1
 * We have to be careful about following facts: 
 *
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where (1-t*t)<<1 is not true
 * Special care has to be taken because of that, which is explained further in the paper. 
 * Here, the CALLING function will have care about always keeping (1-t*t)<<1 .
 *
 * This function is mostly applicable in the low t, low nu regime
 * (Using bessel_integral_lowt_transform should be preferred)
 * */
 void double_bessel_integral_lowt_lownu(double l1 , double l2 , double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Check for large imaginary nu to safely calculate gamma functions
	if(nu_imag>NU_IMAG_LIMIT){
		
		//Calculate total prefactor of 2^(nu-1)
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex_2(1.5 - 0.5*nu_real + (l1-l2)*0.5 , -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2( (l1+l2) + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		
		//Divide result by overflow-safe Gamma(l+3/2)
		double exp_real,exp_imag;
		if((l1+l2)*0.5 < BESSEL_INTEGRAL_MAX_L){
			double scnd_den_gamma = gamma_real( l2 + 1.5);
			exp_factor = exp(lgamma_num_real - a_0+log(t)*l2)/scnd_den_gamma;
			
			exp_real = exp_factor*cos(lgamma_num_imag - b_0);
			exp_imag = exp_factor*sin(lgamma_num_imag - b_0);
		}else{
			double scnd_den_gamma = ln_gamma_real(l2 + 1.5);
			exp_factor = exp(lgamma_num_real - a_0 - scnd_den_gamma+log(t)*l2);
			
			exp_real = exp_factor*cos(lgamma_num_imag - b_0);
			exp_imag = exp_factor*sin(lgamma_num_imag - b_0);
		}
		//We now obtained the total prefactor of the hypergeometric function
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*exp_real - pref_imag*exp_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*exp_real + pref_real*exp_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.5 + (l2-l1)*0.5 - 0.5, nu_imag*0.5, (l1+l2)*0.5 + nu_real*0.5, nu_imag*0.5, l2 + 1.5 ,0.0, t*t, &hypergeom_real, &hypergeom_imag);
		
		//Result
		*res_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
	
	else if ((l1+l2)*0.5 < BESSEL_INTEGRAL_MAX_L){
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate Gamma((3-nu)/2) and Gamma(l+3/2) as the denominator
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real + (l1-l2)*0.5 , -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		double scnd_den_gamma = gamma_real(l2 + 1.5);
		
		//Calculate the numerator
		double gamma_num_real = 0.0; double gamma_num_imag = 0.0;
		gamma_complex( (l1+l2)*0.5 + 0.5*nu_real, 0.5*nu_imag, &gamma_num_real, &gamma_num_imag);
		
		//Calculate fraction of numerator/denominator and multiply with prefactor
		double frac_real = (first_den_gamma_real*gamma_num_real+first_den_gamma_imag*gamma_num_imag)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double frac_imag = (first_den_gamma_real*gamma_num_imag-first_den_gamma_imag*gamma_num_real)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.5 + (l2-l1)*0.5 - 0.5, nu_imag*0.5, (l1+l2)*0.5 + nu_real*0.5, nu_imag*0.5, l2 + 1.5 ,0.0, t*t, &hypergeom_real, &hypergeom_imag);
		
		//Final prefactor of t^l, which is going to make the integral vanishingly small for high l and low t
		double t_l_pref = pow(t,l2);
		*res_real = t_l_pref*(hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = t_l_pref*(hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;		
	}
	else{
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Only Gamma((3-nu+l1-l2)/2) can be calculated directly
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 + (l1-l2)*0.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;

		//Calculate Gamma(l+nu/2)/Gamma(l+3/2)
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2( (l1+l2)*0.5 + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l2 + 1.5);
		exp_factor = exp(lgamma_num_real - scnd_den_gamma); 
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);
		
		//Calculate fraction of numerator and denominator
		double frac_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		double frac_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
		
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.5 + (l2-l1)*0.5 - 0.5, nu_imag*0.5, (l1+l2)*0.5 + nu_real*0.5, nu_imag*0.5, l2 + 1.5 ,0.0, t*t, &hypergeom_real, &hypergeom_imag);
		
		//Final prefactor of t^l, which is going to make the integral vanishingly small for high l and low t
		double t_l_pref = pow(t,l2);
		*res_real = t_l_pref*(hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = t_l_pref*(hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
}

/**
 * This function calculates the integral over the bessel functions given in equation B.1 slightly differently, using a single transformation
 * This transformation reduces the parameters l and nu (divides them by 2), but increases z (roughly factor of 4)
 *
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where (1-t*t)<<1 is not true
 * Special care has to be taken because of that, which is explained further in the paper. 
 * Here, the CALLING function will have care about always keeping (1-t*t)<<1 .
 *
 * This function is mostly applicable in the low t, low nu regime
 * */
void bessel_integral_lowt_transform(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Uses a single transformation to reduce parameters of hypergeo significantly
	//Recurrent parameters
	double z = t*t;
	double z_fac = 4.0*z/(1.0+z)/(1.0+z);
	double log1pz = log(1.0+z);
	//Gamma-overflow safe version
	if(nu_imag>NU_IMAG_LIMIT){
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate Division of Gammas and factor of (1+z)^(l+nu/2)
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex_2(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		exp_factor = exp(lgamma_num_real - a_0 - scnd_den_gamma+log(t)*l-log1pz*(l+0.5*nu_real));
		
		double exp_real = exp_factor*cos(lgamma_num_imag - b_0 - log1pz*0.5*nu_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag - b_0 - log1pz*0.5*nu_imag);
		
		//Multiply this with other prefactor to obtain total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*exp_real - pref_imag*exp_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*exp_real + pref_real*exp_imag);

		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag);
		
		//Result
		*res_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
	else if (l < BESSEL_INTEGRAL_MAX_L){
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate denominator
		//The Gamma((3-nu)/2) and Gamma(l+3/2) can be calculated normally in this case
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		double scnd_den_gamma = gamma_real(l + 1.5);
		
		
		double gamma_num_real = 0.0; double gamma_num_imag = 0.0;
		gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &gamma_num_real, &gamma_num_imag);
		
		
		//Calculate fraction of numerator/denominator and multiply with prefactor
		double frac_real = (first_den_gamma_real*gamma_num_real+first_den_gamma_imag*gamma_num_imag)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double frac_imag = (first_den_gamma_real*gamma_num_imag-first_den_gamma_imag*gamma_num_real)/ first_den_gamma_abs_squared / scnd_den_gamma;
		
		//The old naming convention calls this "tot_pref", but it is not the total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag);
		
		//Calculate the new additional prefactor of (1+z)^(-(l+nu/2))*t^l
		double new_pref_factor = exp(-log1pz*(l+nu_real*0.5) + log(t)*l);
		double new_pref_real = new_pref_factor*cos(-log1pz*nu_imag*0.5);
		double new_pref_imag = new_pref_factor*sin(-log1pz*nu_imag*0.5);
		
		//Multiply old and new prefactor to obtain result
		double bet_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		double bet_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		*res_real = (bet_real*new_pref_real-bet_imag*new_pref_imag);
		*res_imag = (bet_imag*new_pref_real+bet_real*new_pref_imag);
		return;		
	}
	else{
		//Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);

		//Only Gamma((3-nu)/2) can be calculated normally without overflow here
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		
		//Calculate Gamma(l+nu/2)/Gamma(l+3/2)
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		
		exp_factor = exp(lgamma_num_real - scnd_den_gamma);
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);

		//Calculate fraction of numerator and denominator
		double frac_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		double frac_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
		
		//The old naming convention calls this "tot_pref", but it is not the total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);

		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag);

		//Calculate the new additional prefactor of (1+z)^(-(l+nu/2))*t^l
		double new_pref_factor = exp(-log1pz*(l+nu_real*0.5) + log(t)*l);
		double new_pref_real = new_pref_factor*cos(-log1pz*nu_imag*0.5);
		double new_pref_imag = new_pref_factor*sin(-log1pz*nu_imag*0.5);
		//Multiply old and new prefactor to obtain result
		double bet_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		double bet_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		*res_real = (bet_real*new_pref_real-bet_imag*new_pref_imag);
		*res_imag = (bet_imag*new_pref_real+bet_real*new_pref_imag);
		return;	
	}

}
void bessel_integral_lowt_transform_safe(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Uses a single transformation to reduce parameters of hypergeo significantly
	//Recurrent parameters
	double z = t*t;
	double z_fac = 4.0*z/(1.0+z)/(1.0+z);
	double log1pz = log(1.0+z);
	//Gamma-overflow safe version
	if(nu_imag>NU_IMAG_LIMIT){
    int overflows = 0;
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom_safe(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag,&overflows);
		
    //Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate Division of Gammas and factor of (1+z)^(l+nu/2)
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex_2(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		exp_factor = exp(overflows*log(ABS_LIMIT)+lgamma_num_real - a_0 - scnd_den_gamma+log(t)*l-log1pz*(l+0.5*nu_real));
		
		double exp_real = exp_factor*cos(lgamma_num_imag - b_0 - log1pz*0.5*nu_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag - b_0 - log1pz*0.5*nu_imag);
		
		//Multiply this with other prefactor to obtain total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*exp_real - pref_imag*exp_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*exp_real + pref_real*exp_imag);

		
		//Result
		*res_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
	else if (l < BESSEL_INTEGRAL_MAX_L){
    int overflows = 0;
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom_safe(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag, &overflows);
		
    //Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate denominator
		//The Gamma((3-nu)/2) and Gamma(l+3/2) can be calculated normally in this case
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		double scnd_den_gamma = gamma_real(l + 1.5);
		
		
		double gamma_num_real = 0.0; double gamma_num_imag = 0.0;
		gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &gamma_num_real, &gamma_num_imag);
		
		
		//Calculate fraction of numerator/denominator and multiply with prefactor
		double frac_real = (first_den_gamma_real*gamma_num_real+first_den_gamma_imag*gamma_num_imag)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double frac_imag = (first_den_gamma_real*gamma_num_imag-first_den_gamma_imag*gamma_num_real)/ first_den_gamma_abs_squared / scnd_den_gamma;
		
		//The old naming convention calls this "tot_pref", but it is not the total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		
		//Calculate the new additional prefactor of (1+z)^(-(l+nu/2))*t^l
		double new_pref_factor = exp(-log1pz*(l+nu_real*0.5) + log(t)*l + overflows*log(ABS_LIMIT));
		double new_pref_real = new_pref_factor*cos(-log1pz*nu_imag*0.5);
		double new_pref_imag = new_pref_factor*sin(-log1pz*nu_imag*0.5);
		
		//Multiply old and new prefactor to obtain result
		double bet_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		double bet_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		*res_real = (bet_real*new_pref_real-bet_imag*new_pref_imag);
		*res_imag = (bet_imag*new_pref_real+bet_real*new_pref_imag);
		return;		
	}
	else{
    int overflows = 0;
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom_safe(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag,&overflows);
		
    //Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);

		//Only Gamma((3-nu)/2) can be calculated normally without overflow here
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		
		//Calculate Gamma(l+nu/2)/Gamma(l+3/2)
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		
		exp_factor = exp(lgamma_num_real - scnd_den_gamma);
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);

		//Calculate fraction of numerator and denominator
		double frac_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		double frac_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
		
		//The old naming convention calls this "tot_pref", but it is not the total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);

		//Calculate the new additional prefactor of (1+z)^(-(l+nu/2))*t^l
		double new_pref_factor = exp(-log1pz*(l+nu_real*0.5) + log(t)*l + overflows*log(ABS_LIMIT));
		double new_pref_real = new_pref_factor*cos(-log1pz*nu_imag*0.5);
		double new_pref_imag = new_pref_factor*sin(-log1pz*nu_imag*0.5);
		
		//Multiply old and new prefactor to obtain result
		double bet_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		double bet_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		*res_real = (bet_real*new_pref_real-bet_imag*new_pref_imag);
		*res_imag = (bet_imag*new_pref_real+bet_real*new_pref_imag);
		return;	
	}

}
void bessel_integral_lowt_transform_safe_flag(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag,short* overflow_flag){
	//Uses a single transformation to reduce parameters of hypergeo significantly
	//Recurrent parameters
	double z = t*t;
	double z_fac = 4.0*z/(1.0+z)/(1.0+z);
	double log1pz = log(1.0+z);
	//Gamma-overflow safe version
	if(nu_imag>NU_IMAG_LIMIT){
    int overflows = 0;
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom_safe(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag,&overflows);
		if(overflows>0){*overflow_flag = _TRUE_;}else{*overflow_flag=_FALSE_;}
    //Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate Division of Gammas and factor of (1+z)^(l+nu/2)
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex_2(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		exp_factor = exp(overflows*log(ABS_LIMIT)+lgamma_num_real - a_0 - scnd_den_gamma+log(t)*l-log1pz*(l+0.5*nu_real));
		
		double exp_real = exp_factor*cos(lgamma_num_imag - b_0 - log1pz*0.5*nu_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag - b_0 - log1pz*0.5*nu_imag);
		
		//Multiply this with other prefactor to obtain total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*exp_real - pref_imag*exp_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*exp_real + pref_real*exp_imag);

		
		//Result
		*res_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		*res_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		return;
	}
	else if (l < BESSEL_INTEGRAL_MAX_L){
    int overflows = 0;
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom_safe(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag, &overflows);
		if(overflows>0){*overflow_flag = _TRUE_;}else{*overflow_flag=_FALSE_;}
		
    //Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);
		
		//Calculate denominator
		//The Gamma((3-nu)/2) and Gamma(l+3/2) can be calculated normally in this case
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		double scnd_den_gamma = gamma_real(l + 1.5);
		
		
		double gamma_num_real = 0.0; double gamma_num_imag = 0.0;
		gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &gamma_num_real, &gamma_num_imag);
		
		
		//Calculate fraction of numerator/denominator and multiply with prefactor
		double frac_real = (first_den_gamma_real*gamma_num_real+first_den_gamma_imag*gamma_num_imag)/ first_den_gamma_abs_squared / scnd_den_gamma;
		double frac_imag = (first_den_gamma_real*gamma_num_imag-first_den_gamma_imag*gamma_num_real)/ first_den_gamma_abs_squared / scnd_den_gamma;
		
		//The old naming convention calls this "tot_pref", but it is not the total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);
		
		
		//Calculate the new additional prefactor of (1+z)^(-(l+nu/2))*t^l
		double new_pref_factor = exp(-log1pz*(l+nu_real*0.5) + log(t)*l + overflows*log(ABS_LIMIT));
		double new_pref_real = new_pref_factor*cos(-log1pz*nu_imag*0.5);
		double new_pref_imag = new_pref_factor*sin(-log1pz*nu_imag*0.5);
		
		//Multiply old and new prefactor to obtain result
		double bet_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		double bet_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		*res_real = (bet_real*new_pref_real-bet_imag*new_pref_imag);
		*res_imag = (bet_imag*new_pref_real+bet_real*new_pref_imag);
		return;		
	}
	else{
    int overflows = 0;
		//Calculate hypergeometric function
		double hypergeom_real = 0.0; double hypergeom_imag = 0.0;
		hypergeom_safe(nu_real*0.25 + 0.5*l, nu_imag*0.25, l*0.5+0.5 + nu_real*0.25, nu_imag*0.25, l + 1.5 ,0.0, z_fac, &hypergeom_real, &hypergeom_imag,&overflows);
		if(overflows>0){*overflow_flag = _TRUE_;}else{*overflow_flag=_FALSE_;}
		
    //Calculate total prefactor of hypergeometric function
		double exp_factor = exp(LOG_2*nu_real - LOG_2);
		double pref_real = exp_factor*cos(LOG_2*nu_imag);
		double pref_imag = exp_factor*sin(LOG_2*nu_imag);

		//Only Gamma((3-nu)/2) can be calculated normally without overflow here
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		
		//Calculate Gamma(l+nu/2)/Gamma(l+3/2)
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		
		exp_factor = exp(lgamma_num_real - scnd_den_gamma);
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);

		//Calculate fraction of numerator and denominator
		double frac_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		double frac_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
		
		//The old naming convention calls this "tot_pref", but it is not the total prefactor
		double tot_pref_real = MATH_PI*MATH_PI* (pref_real*frac_real - pref_imag*frac_imag);
		double tot_pref_imag = MATH_PI*MATH_PI* (pref_imag*frac_real + pref_real*frac_imag);

		//Calculate the new additional prefactor of (1+z)^(-(l+nu/2))*t^l
		double new_pref_factor = exp(-log1pz*(l+nu_real*0.5) + log(t)*l + overflows*log(ABS_LIMIT));
		double new_pref_real = new_pref_factor*cos(-log1pz*nu_imag*0.5);
		double new_pref_imag = new_pref_factor*sin(-log1pz*nu_imag*0.5);
		
		//Multiply old and new prefactor to obtain result
		double bet_real = (hypergeom_real*tot_pref_real - hypergeom_imag*tot_pref_imag);
		double bet_imag = (hypergeom_real*tot_pref_imag + hypergeom_imag*tot_pref_real);
		*res_real = (bet_real*new_pref_real-bet_imag*new_pref_imag);
		*res_imag = (bet_imag*new_pref_real+bet_real*new_pref_imag);
		return;	
	}

}
/**
 * This function calculates the integral over the bessel functions given in equation B.6
 * 
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1
 * Special care has to be taken because of that, which is explained further in the paper. 
 * Here, the CALLING function will have care about always keeping t<<1 .
 *
 * This function is mostly applicable in the high t, low nu regime
 * (Using bessel_integral_hight_transform should be preferred)
 * */
void bessel_integral_hight_lownu(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Gamma-overflow safe evaluations
	if(nu_imag>NU_IMAG_LIMIT){
		//Needed parameter t_tilde
		double t_tilde = -(1 - t*t)*(1 - t*t) / (4 * t*t);
		
		//Prefactor of t^(-nu/2)*pi^(3/2)
		double exp_factor = pow(t, -nu_real*0.5);
		double temp = -log(t)*nu_imag*0.5;
		double tot_pref_real = exp_factor*cos(temp)*MATH_PI_3_2;
		double tot_pref_imag = exp_factor*sin(temp)*MATH_PI_3_2;
		
		//Calculate loggammas
		double a_0 = 0.0; double b_0 = 0.0;
		ln_gamma_complex_2(1.5 - nu_real*0.5, -nu_imag*0.5, &a_0, &b_0);
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		ln_gamma_complex_2(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
		double a_4 = 0.0; double b_4 = 0.0;
		ln_gamma_complex_2(nu_real*0.5 - 1, nu_imag*0.5, &a_4, &b_4);
		
		//Calculate prefactor of first hypergeom
		exp_factor = exp(a_1 + a_2 - a_3 - a_0);
		double pref_first_real = exp_factor*cos(b_1 + b_2 - b_3 - b_0);
		double pref_first_imag = exp_factor*sin(b_1 + b_2 - b_3 - b_0);
		
		//Calculate part of prefactor of second hypergeom
		exp_factor = exp(a_4 - a_0);
		double pref_scnd_gamma_real = exp_factor*cos(b_4 - b_0);
		double pref_scnd_gamma_imag = exp_factor*sin(b_4 - b_0);
		
		//Calculate first hypergeom
		double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
		hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 - l*0.5 - 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, t_tilde, &first_hypergeom_real, &first_hypergeom_imag);
		
		//Calculate additional prefacotr of second hypergeom
		exp_factor = pow(-t_tilde*0.25, 1 - nu_real*0.5);
		temp = -log(-t_tilde*0.25)*nu_imag*0.5;
		
		double pref_scnd_real = exp_factor*cos(temp);
		double pref_scnd_imag = exp_factor*sin(temp);
		
		//Calculate final prefactor of second hypergeom
		double tot_pref_scnd_real = pref_scnd_real*pref_scnd_gamma_real - pref_scnd_imag*pref_scnd_gamma_imag;
		double tot_pref_scnd_imag = pref_scnd_real*pref_scnd_gamma_imag + pref_scnd_imag*pref_scnd_gamma_real;
		
		//Calculate second hypergeom
		double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
		hypergeom(l*0.5 - nu_real*0.25 + 1, -nu_imag*0.25, 0.5 - l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, t_tilde, &scnd_hypergeom_real, &scnd_hypergeom_imag);
		
		//Calculate sum of hypergeoms with their prefactors
		double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
		double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
		
		//Calculate total
		*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
		*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
		return;	
	}
	else if (l<BESSEL_INTEGRAL_MAX_L){
		//Calculate t_tilde
		double t_tilde = -(1 - t*t)*(1 - t*t) / (4 * t*t);
		
		//Calculate prefactor of t^(-nu/2)*pi^(3/2)
		double exp_factor = pow(t, -nu_real*0.5);
		double temp = -log(t)*nu_imag*0.5;
		double pref_real = exp_factor*cos(temp)*MATH_PI_3_2;
		double pref_imag = exp_factor*sin(temp)*MATH_PI_3_2;
		
		double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
		gamma_complex(1.5 - nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
		
		double pref_den_gamma_abs_squared = pref_den_gamma_real*pref_den_gamma_real + pref_den_gamma_imag*pref_den_gamma_imag;
		double tot_pref_real = (pref_real*pref_den_gamma_real + pref_imag*pref_den_gamma_imag) / pref_den_gamma_abs_squared;
		double tot_pref_imag = (pref_imag*pref_den_gamma_real - pref_real*pref_den_gamma_imag) / pref_den_gamma_abs_squared;
		
		//Calculate prefactor of first hypergeom
		double a_1 = 0.0; double b_1 = 0.0;
		gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
		double abs_3 = (a_3*a_3 + b_3*b_3);
		
		double pref_first_real = (a_1*a_2*a_3 - b_1*b_2*a_3 + b_1*b_3*a_2 + b_2*b_3*a_1) / abs_3;
		double pref_first_imag = (b_1*a_2*a_3 + b_2*a_1*a_3 - b_3*a_1*a_2 + b_1*b_2*b_3) / abs_3;
		
		//Calculate first hypergeom
		double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
		hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 - l*0.5 - 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, t_tilde, &first_hypergeom_real, &first_hypergeom_imag);
		
		//Calculate prefactor of second hypergeom
		exp_factor = pow(-t_tilde*0.25, 1 - nu_real*0.5);
		temp = -log(-t_tilde*0.25)*nu_imag*0.5;
		double pref_scnd_real = exp_factor*cos(temp);
		double pref_scnd_imag = exp_factor*sin(temp);
		
		double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
		gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
		
		double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
		double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
		
		//Calculate second hypergeom
		double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
		hypergeom(l*0.5 - nu_real*0.25 + 1, -nu_imag*0.25, 0.5 - l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, t_tilde, &scnd_hypergeom_real, &scnd_hypergeom_imag);
		
		
		
		//Calculate sum of hypergeoms with their prefactors
		double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
		double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
		
		//Calculate total
		*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
		*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
		return;
	}
	else{
		//Calculate t_tilde
		double t_tilde = -(1 - t*t)*(1 - t*t) / (4 * t*t);
		
		//Calculate prefactor of t^(-nu/2)*pi^(3/2)
		double exp_factor = pow(t, -nu_real*0.5);
		double temp = -log(t)*nu_imag*0.5;
		double pref_real = exp_factor*cos(temp)*MATH_PI_3_2;
		double pref_imag = exp_factor*sin(temp)*MATH_PI_3_2;
		
		//Calculate first prefactor
		double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
		gamma_complex(1.5 - nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
		
		double pref_den_gamma_abs_squared = pref_den_gamma_real*pref_den_gamma_real + pref_den_gamma_imag*pref_den_gamma_imag;
		double tot_pref_real = (pref_real*pref_den_gamma_real + pref_imag*pref_den_gamma_imag) / pref_den_gamma_abs_squared;
		double tot_pref_imag = (pref_imag*pref_den_gamma_real - pref_real*pref_den_gamma_imag) / pref_den_gamma_abs_squared;
		
		//Calculate prefactor of first hypergeom
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
		exp_factor = exp(a_1 - a_3);
		double exp_real = exp_factor*cos(b_1 - b_3);
		double exp_imag = exp_factor*sin(b_1 - b_3);
		
		double pref_first_real = a_2*exp_real - b_2*exp_imag;
		double pref_first_imag = a_2*exp_imag + b_2*exp_real;
		
		//Calculate first hypergeom
		double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
		hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 - l*0.5 - 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, t_tilde, &first_hypergeom_real, &first_hypergeom_imag);
		
		//Calculate prefactor of second hypergeom
		exp_factor = pow(-t_tilde*0.25, 1 - nu_real*0.5);
		temp = -log(-t_tilde*0.25)*nu_imag*0.5;
		double pref_scnd_real = exp_factor*cos(temp);
		double pref_scnd_imag = exp_factor*sin(temp);
		
		double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
		gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
		double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
		double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
		
		//Calculate second hypergeom
		double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
		hypergeom(l*0.5 - nu_real*0.25 + 1, -nu_imag*0.25, 0.5 - l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, t_tilde, &scnd_hypergeom_real, &scnd_hypergeom_imag);
		
		//Calculate sum of hypergeoms with their prefactors
		double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
		double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
		
		//Calculate total
		*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
		*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
		return;
	}
}
/**
 * This function calculates the integral over the bessel functions given in equation B.6, but with two different l
 * 
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1
 * Special care has to be taken because of that, which is explained further in the paper. 
 * Here, the CALLING function will have care about always keeping t<<1 .
 *
 * This function is mostly applicable in the high t, low nu regime
 * (Using bessel_integral_hight_transform should be preferred)
 * */
void double_bessel_integral_hight_lownu(double l1, double l2, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
  double first_pref_real; double first_pref_imag;
  double second_pref_real; double second_pref_imag;
  if( nu_imag > NU_IMAG_LIMIT){
    double a_0 = 0.0; double b_0 = 0.0;
		ln_gamma_complex_2(nu_real - 2, nu_imag, &a_0, &b_0);
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex_2(( 3 - nu_real + l1 - l2)*0.5, -nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		ln_gamma_complex_2((nu_real - 1 - l1 + l2)*0.5, nu_imag*0.5, &a_2, &b_2);
		
    double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1)+a_0-a_1-a_2);
    first_pref_real = exp_factor*cos(LOG_2*nu_imag+b_0-b_1-b_2);
    first_pref_imag = exp_factor*sin(LOG_2*nu_imag+b_0-b_1-b_2);
  
  }
  else{
    double num_gamma_real = 0.0; double num_gamma_imag = 0.0;
    gamma_complex( nu_real - 2, nu_imag, &num_gamma_real, &num_gamma_imag);
    double den1_gamma_real = 0.0; double den1_gamma_imag = 0.0;
    gamma_complex(( 3 - nu_real + l1 - l2)*0.5, -nu_imag*0.5, &den1_gamma_real, &den1_gamma_imag);
    double den2_gamma_real = 0.0; double den2_gamma_imag = 0.0;
    gamma_complex((nu_real - 1 - l1 + l2)*0.5, nu_imag*0.5, &den2_gamma_real, &den2_gamma_imag);
    
    double den_real = den1_gamma_real*den2_gamma_real - den1_gamma_imag*den2_gamma_imag;
    double den_imag = den1_gamma_real*den2_gamma_imag + den1_gamma_imag*den2_gamma_real;
    
    double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1));
    double exp_real = exp_factor*cos(LOG_2*nu_imag);
    double exp_imag = exp_factor*sin(LOG_2*nu_imag);
    double num_real = num_gamma_real*exp_real - num_gamma_imag*exp_imag;
    double num_imag = num_gamma_real*exp_imag + num_gamma_imag*exp_real;
    
    double den_abs_squared = den_real*den_real + den_imag*den_imag;
    first_pref_real = (num_real*den_real+num_imag*den_imag)/den_abs_squared;
    first_pref_imag = (num_imag*den_real-num_real*den_imag)/den_abs_squared;
    
  }
  //Prefactor 1 is calculated
  if(nu_imag > NU_IMAG_LIMIT){
    double a_0 = 0.0; double b_0 = 0.0;
		ln_gamma_complex_2((l1 + l2 + nu_real)*0.5, nu_imag*0.5, &a_0, &b_0);
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex_2(2 - nu_real, nu_imag, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		ln_gamma_complex_2((l2+l1 + 4 - nu_real)*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		ln_gamma_complex_2((l2-l1 + 3 - nu_real)*0.5, -nu_imag*0.5, &a_3, &b_3);
		double a_4 = 0.0; double b_4 = 0.0;
		ln_gamma_complex_2((l1-l2 + 3 - nu_real)*0.5, -nu_imag*0.5, &a_4, &b_4);
		
    double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1)+a_0+a_1-a_2-a_3-a_4);
    second_pref_real = exp_factor*cos(LOG_2*nu_imag+b_0+b_1-b_2-b_3-b_4);
    second_pref_imag = exp_factor*sin(LOG_2*nu_imag+b_0+b_1-b_2-b_3-b_4);
    
  }
  else{
    double bet_real; double bet_imag;
    if((l1+l2)*0.5 > BESSEL_INTEGRAL_MAX_L ) {
      double a_0 = 0.0; double b_0 = 0.0;
      ln_gamma_complex_2((l1 + l2 + nu_real)*0.5, nu_imag*0.5, &a_0, &b_0);
      double a_2 = 0.0; double b_2 = 0.0;
      ln_gamma_complex_2((l2+l1 + 4 - nu_real)*0.5, -nu_imag*0.5, &a_2, &b_2);
      double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1)+a_0-a_2);
      bet_real = exp_factor*cos(LOG_2*nu_imag+b_0-b_2);
      bet_imag = exp_factor*sin(LOG_2*nu_imag+b_0-b_2);
    
    }else{
      double n_gamma_real = 0.0; double n_gamma_imag = 0.0;
      gamma_complex((l1 + l2 + nu_real)*0.5,nu_imag*0.5,&n_gamma_real,&n_gamma_imag);
      double d_gamma_real = 0.0; double d_gamma_imag = 0.0;
      gamma_complex((l2 + l1 +4 -nu_real)*0.5,-nu_imag*0.5,&d_gamma_real,&d_gamma_imag);
      
      double d_gamma_abs_squared = d_gamma_real*d_gamma_real + d_gamma_imag*d_gamma_imag;
      double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1));
      double exp_real = exp_factor*cos(LOG_2*nu_imag);
      double exp_imag = exp_factor*sin(LOG_2*nu_imag);
      
      double num_real = n_gamma_real*exp_real - n_gamma_imag*exp_imag;
      double num_imag = n_gamma_real*exp_imag + n_gamma_imag*exp_real;
      
      bet_real = ( num_real*d_gamma_real + num_imag*d_gamma_imag)/d_gamma_abs_squared;
      bet_imag = ( num_imag*d_gamma_real - num_real*d_gamma_imag)/d_gamma_abs_squared;
      
    }
    
    
    
    double num_gamma_real = 0.0; double num_gamma_imag = 0.0;
    gamma_complex(2-nu_real,-nu_imag,&num_gamma_real,&num_gamma_imag);
    double den1_gamma_real = 0.0; double den1_gamma_imag = 0.0;
    gamma_complex((l2+3-l1-nu_real)*0.5,-nu_imag*0.5,&den1_gamma_real,&den1_gamma_imag);
    double den2_gamma_real = 0.0; double den2_gamma_imag = 0.0;
    gamma_complex((l2+3-l1-nu_real)*0.5,-nu_imag*0.5,&den2_gamma_real,&den2_gamma_imag);
    
    double den_real = den1_gamma_real*den2_gamma_real - den1_gamma_imag*den2_gamma_imag;
    double den_imag = den1_gamma_real*den2_gamma_imag + den1_gamma_imag*den2_gamma_real;
    
    double num_real = num_gamma_real*bet_real - num_gamma_imag*bet_imag;
    double num_imag = num_gamma_real*bet_imag + num_gamma_imag*bet_real;
    
    second_pref_real = num_real*den_real+num_imag*den_imag;
    second_pref_imag = num_imag*den_real-num_real*den_imag;
  }
  //Needed parameter t_tilde
  double t_tilde = (1 - t*t)*(1 - t*t);
  
  double exp_factor = pow(t_tilde, 2-nu_real);
  double temp = -log(t_tilde)*nu_imag;
  double tot_pref_real = exp_factor*cos(temp);
  double tot_pref_imag = exp_factor*sin(temp);
    
  double tot_pref_first_real = tot_pref_real*first_pref_real - tot_pref_imag*first_pref_imag;
  double tot_pref_first_imag = tot_pref_imag*first_pref_real + tot_pref_real*first_pref_imag;
  
  //Calculate first hypergeom
  double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
  hypergeom((l2+3-l1-nu_real)*0.5, -nu_imag*0.5, (l2+l1+3-nu_real)*0.5, -nu_imag*0.5, 3-nu_real, -nu_imag, t_tilde, &first_hypergeom_real, &first_hypergeom_imag);
  
  //Calculate second hypergeom
  double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
  hypergeom((l2-1-l1+nu_real)*0.5, nu_imag*0.5, (l1+l2+nu_real)*0.5, nu_imag*0.5, nu_real-1, nu_imag, t_tilde, &scnd_hypergeom_real, &scnd_hypergeom_imag);
  
  //Calculate sum of hypergeoms with their prefactors
  double total_bracket_real = (first_hypergeom_real*second_pref_real - first_hypergeom_imag*second_pref_imag) + (scnd_hypergeom_real*tot_pref_first_real - scnd_hypergeom_imag*tot_pref_first_imag);
  double total_bracket_imag = (first_hypergeom_real*second_pref_imag + first_hypergeom_imag*second_pref_real) + (scnd_hypergeom_real*tot_pref_first_imag + scnd_hypergeom_imag*tot_pref_first_real);
  
  //Calculate total
  *res_real = total_bracket_real;
  *res_imag = total_bracket_imag;
  return;	
}

/**
 * This function calculates the integral over the bessel functions given in equation B.6 in a transformed version 
 * ( with a different derivation, but similar form )
 * 
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1
 * Special care has to be taken because of that, which is explained further in the paper. 
 * Here, the CALLING function will have care about always keeping t<<1 .
 *
 * This function is mostly applicable in the high t regime
 * */
void bessel_integral_hight_transform(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Recurrent transformed factors of z
	double z = t*t;
  //printf("High t = %.10e ||| ",t);
	double z_fac1 = ((1-z)*(1-z))/((1+z)*(1+z));
	double z_fac2 = 2*(1+z)/(1-z);
	double lnz_fac3 = log(2.0/(1.0+z));
  //printf("Calculating with argument %.10e \n",z_fac1);
	if(nu_imag<NU_IMAG_LIMIT){
		if(l<BESSEL_INTEGRAL_MAX_L){
			double exp_factor = exp(lnz_fac3*(l+nu_real*0.5));
			double exp_phase = lnz_fac3*nu_imag*0.5;
			double exp_real = exp_factor*cos(exp_phase);
			double exp_imag = exp_factor*sin(exp_phase);
			
			//Calculate total prefactor
			double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
			gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
			double pref_den_gamma_abs_squared = (pref_den_gamma_real*pref_den_gamma_real+pref_den_gamma_imag*pref_den_gamma_imag);
			double tot_pref_real = pow(t,l)*MATH_PI_3_2*(exp_real*pref_den_gamma_real+exp_imag*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
			double tot_pref_imag = pow(t,l)*MATH_PI_3_2*(exp_imag*pref_den_gamma_real-exp_real*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
		
			double a_1 = 0.0; double b_1 = 0.0;
			gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
			double a_2 = 0.0; double b_2 = 0.0;
			gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
			double a_3 = 0.0; double b_3 = 0.0;
			gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
			double abs_3 = a_3*a_3 + b_3*b_3;
			exp_real = (a_1*a_3+b_1*b_3)/abs_3;
			exp_imag = (b_1*a_3-a_1*b_3)/abs_3;
			double pref_first_real = a_2*exp_real - b_2*exp_imag;
			double pref_first_imag = a_2*exp_imag + b_2*exp_real;
			
			//Calculate first hypergeom
			double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
			hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag);
			
			//Calculate prefactor of second hypergeom
			exp_factor = pow(z_fac2, nu_real-2.0);
			double temp = log(z_fac2)*nu_imag;
			double pref_scnd_real = exp_factor*cos(temp);
			double pref_scnd_imag = exp_factor*sin(temp);
			
			double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
			gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
			double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
			double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
			
			//Calculate second hypergeom
			double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
			hypergeom(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag);
      
      //Calculate sum of hypergeoms with their prefactors
			double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
			double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
			
			//Calculate total
			*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
			*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
			return;
		}
		else{
			//Calculate prefactor including t^l
			double exp_factor = exp(lnz_fac3*(l+nu_real*0.5)+log(t)*l);
			double exp_phase = lnz_fac3*nu_imag*0.5;
			double exp_real = exp_factor*cos(exp_phase);
			double exp_imag = exp_factor*sin(exp_phase);
			
			//Calculate total prefactor
			double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
			gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
			double pref_den_gamma_abs_squared = (pref_den_gamma_real*pref_den_gamma_real+pref_den_gamma_imag*pref_den_gamma_imag);
			double tot_pref_real = MATH_PI_3_2*(exp_real*pref_den_gamma_real+exp_imag*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
			double tot_pref_imag = MATH_PI_3_2*(exp_imag*pref_den_gamma_real-exp_real*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
			//printf("Tot pref = %.10e+%.10ej \n",tot_pref_real,tot_pref_imag);
      
			//Calculate gamma factors for large l 
			double a_1 = 0.0; double b_1 = 0.0;
			ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
			double a_2 = 0.0; double b_2 = 0.0;
			gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
			double a_3 = 0.0; double b_3 = 0.0;
			ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
			
			exp_factor = exp(a_1 - a_3);
			exp_real = exp_factor*cos(b_1 - b_3);
			exp_imag = exp_factor*sin(b_1 - b_3);
			//printf("Gamma prod = %.10e+%.10ej \n",exp_real,exp_imag);
      
			//Calculate first prefactor
			double pref_first_real = a_2*exp_real - b_2*exp_imag;
			double pref_first_imag = a_2*exp_imag + b_2*exp_real;
			//printf("Pref first = %.15e+%.15ej \n",pref_first_real,pref_first_imag);
			//Calculate first hypergeom
			double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
			hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag);
			
			
			//Calculate prefactor of second hypergeom
			exp_factor = pow(z_fac2, nu_real-2.0);
			double temp = log(z_fac2)*nu_imag;
			double pref_scnd_real = exp_factor*cos(temp);
			double pref_scnd_imag = exp_factor*sin(temp);
			
			double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
			gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
			double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
			double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
			
			//Calculate second hypergeom
			double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
			hypergeom(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag);
			
			//printf("Pref Scnd = %.15e+%.15ej \n",tot_pref_scnd_real,tot_pref_scnd_imag);
			//Calculate sum of hypergeoms with their prefactors
			double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
			double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
			//printf("Bracket = %.10e+%.10ej \n",total_bracket_real,total_bracket_imag);
			//Calculate total
			*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
			*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
			return;
		}
	}else{
		//Calculate prefactor including t^l
		double exp_factor = exp(lnz_fac3*(l+nu_real*0.5)+log(t)*l);
		double exp_phase = lnz_fac3*nu_imag*0.5;
		double exp_real = exp_factor*cos(exp_phase);
		double exp_imag = exp_factor*sin(exp_phase);
		
		double tot_pref_real = MATH_PI_3_2*exp_real;
		double tot_pref_imag = MATH_PI_3_2*exp_imag;
		
		//Calculate gamma factors (large nu -> lnGamma)
		double a_0 = 0.0; double b_0 = 0.0;
		ln_gamma_complex_2(1.5-nu_real*0.5, -nu_imag*0.5, &a_0, &b_0);
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		ln_gamma_complex_2(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);

		exp_factor = exp(a_1+a_2 - a_3 -a_0);
		exp_real = exp_factor*cos(b_1+b_2 - b_3-b_0);
		exp_imag = exp_factor*sin(b_1+b_2 - b_3-b_0);
		
		//Calculate first prefactor
		double pref_first_real = exp_real;
		double pref_first_imag = exp_imag;
		
		//Calculate first hypergeom
		double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
		hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag);
		
		//Calculate prefactor of second hypergeom
		exp_factor = pow(z_fac2, nu_real-2.0);
		double temp = log(z_fac2)*nu_imag;
		double pref_scnd_real = exp_factor*cos(temp);
		double pref_scnd_imag = exp_factor*sin(temp);
		
		double c_0 = 0.0; double d_0 = 0.0;
		ln_gamma_complex_2(nu_real*0.5 - 1, nu_imag*0.5, &c_0, &d_0);
		exp_factor = exp(c_0-a_0);
		exp_real = exp_factor*cos(d_0-b_0);
		exp_imag = exp_factor*sin(d_0-b_0);
		
		//Calculate second prefactor
		double pref_num_gamma_real = exp_real;
		double pref_num_gamma_imag = exp_imag;
		double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
		double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
		
		//Calculate second hypergeom
		double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
		hypergeom(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag);
		
		//Calculate sum of hypergeoms with their prefactors
		double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
		double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
		
		//Calculate total
		*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
		*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
		return;	
	}
}
/*
void bessel_integral_hight_transform_safe(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Recurrent transformed factors of z
	double z = t*t;
  //printf("High t = %.10e ||| ",t);
	double z_fac1 = ((1-z)*(1-z))/((1+z)*(1+z));
	double z_fac2 = 2*(1+z)/(1-z);
	double lnz_fac3 = log(2.0/(1.0+z));
  //printf("Calculating with argument %.10e \n",z_fac1);
	if(nu_imag<NU_IMAG_LIMIT){
		if(l<BESSEL_INTEGRAL_MAX_L){
			int overflows_1 = 0; 
			//Calculate first hypergeom
			double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
			hypergeom_safe(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag,&overflows_1);
			
      int overflows_2 = 0;
      //Calculate second hypergeom
			double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
			hypergeom_safe(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag,&overflows_2);
			
      while(overflows_1 > overflows_2){
        first_hypergeom_real/=ABS_LIMIT;
        first_hypergeom_imag/=ABS_LIMIT;
        //The two functions are very close to each other, so this should never happen
      }
      while(overflows_2 > overflows_1){
        scnd_hypergeom_real/=ABS_LIMIT;
        scnd_hypergeom_imag/=ABS_LIMIT;
        //The two functions are very close to each other, so this should never happen
      }
      double exp_factor = exp(lnz_fac3*(l+nu_real*0.5)+log(t)*l-overflows_1*log(ABS_LIMIT));
			double exp_phase = lnz_fac3*nu_imag*0.5;
			double exp_real = exp_factor*cos(exp_phase);
			double exp_imag = exp_factor*sin(exp_phase);
			
			//Calculate total prefactor
			double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
			gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
			double pref_den_gamma_abs_squared = (pref_den_gamma_real*pref_den_gamma_real+pref_den_gamma_imag*pref_den_gamma_imag);
			double tot_pref_real = MATH_PI_3_2*(exp_real*pref_den_gamma_real+exp_imag*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
			double tot_pref_imag = MATH_PI_3_2*(exp_imag*pref_den_gamma_real-exp_real*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
		
			double a_1 = 0.0; double b_1 = 0.0;
			gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
			double a_2 = 0.0; double b_2 = 0.0;
			gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
			double a_3 = 0.0; double b_3 = 0.0;
			gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
			double abs_3 = a_3*a_3 + b_3*b_3;
			exp_real = (a_1*a_3+b_1*b_3)/abs_3;
			exp_imag = (b_1*a_3-a_1*b_3)/abs_3;
			double pref_first_real = a_2*exp_real - b_2*exp_imag;
			double pref_first_imag = a_2*exp_imag + b_2*exp_real;
			
			//Calculate prefactor of second hypergeom
			exp_factor = pow(z_fac2, nu_real-2.0);
			double temp = log(z_fac2)*nu_imag;
			double pref_scnd_real = exp_factor*cos(temp);
			double pref_scnd_imag = exp_factor*sin(temp);
			
			double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
			gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
			double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
			double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
			
			//Calculate sum of hypergeoms with their prefactors
			double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
			double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
			
			//Calculate total
			*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
			*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
			return;
		}
		else{
      int overflows_1 = 0;
      //Calculate first hypergeom
			double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
			hypergeom_safe(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag,&overflows_1);
			
      int overflows_2 = 0;
			//Calculate second hypergeom
			double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
			hypergeom_safe(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag,&overflows_2);
			
      while(overflows_1 > overflows_2){
        first_hypergeom_real/=ABS_LIMIT;
        first_hypergeom_imag/=ABS_LIMIT;
        //The two functions are very close to each other, so this should never happen
      }
      while(overflows_2 > overflows_1){
        scnd_hypergeom_real/=ABS_LIMIT;
        scnd_hypergeom_imag/=ABS_LIMIT;
        //The two functions are very close to each other, so this should never happen
      }
      
			//Calculate prefactor including t^l
			double exp_factor = exp(lnz_fac3*(l+nu_real*0.5)+log(t)*l-overflows_1*log(ABS_LIMIT));
			double exp_phase = lnz_fac3*nu_imag*0.5;
			double exp_real = exp_factor*cos(exp_phase);
			double exp_imag = exp_factor*sin(exp_phase);
			
			//Calculate total prefactor
			double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
			gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
			double pref_den_gamma_abs_squared = (pref_den_gamma_real*pref_den_gamma_real+pref_den_gamma_imag*pref_den_gamma_imag);
			double tot_pref_real = MATH_PI_3_2*(exp_real*pref_den_gamma_real+exp_imag*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
			double tot_pref_imag = MATH_PI_3_2*(exp_imag*pref_den_gamma_real-exp_real*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
			//printf("Tot pref = %.10e+%.10ej \n",tot_pref_real,tot_pref_imag);
      
			//Calculate gamma factors for large l 
			double a_1 = 0.0; double b_1 = 0.0;
			ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
			double a_2 = 0.0; double b_2 = 0.0;
			gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
			double a_3 = 0.0; double b_3 = 0.0;
			ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
			
			exp_factor = exp(a_1 - a_3);
			exp_real = exp_factor*cos(b_1 - b_3);
			exp_imag = exp_factor*sin(b_1 - b_3);
			//printf("Gamma prod = %.10e+%.10ej \n",exp_real,exp_imag);
      
			//Calculate first prefactor
			double pref_first_real = a_2*exp_real - b_2*exp_imag;
			double pref_first_imag = a_2*exp_imag + b_2*exp_real;
			//printf("Pref first = %.15e+%.15ej \n",pref_first_real,pref_first_imag);
			
			//Calculate prefactor of second hypergeom
			exp_factor = pow(z_fac2, nu_real-2.0);
			double temp = log(z_fac2)*nu_imag;
			double pref_scnd_real = exp_factor*cos(temp);
			double pref_scnd_imag = exp_factor*sin(temp);
			
			double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
			gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
			double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
			double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
			
			//printf("Pref Scnd = %.15e+%.15ej \n",tot_pref_scnd_real,tot_pref_scnd_imag);
			//Calculate sum of hypergeoms with their prefactors
			double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
			double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
			//printf("Bracket = %.10e+%.10ej \n",total_bracket_real,total_bracket_imag);
			//Calculate total
			*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
			*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
			return;
		}
	}else{
    int overflows_1 = 0;
    //Calculate first hypergeom
		double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
		hypergeom_safe(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag,&overflows_1);
		
    int overflows_2 = 0;
		//Calculate second hypergeom
		double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
		hypergeom_safe(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag, &overflows_2);
		
    while(overflows_1 > overflows_2){
      first_hypergeom_real/=ABS_LIMIT;
      first_hypergeom_imag/=ABS_LIMIT;
      //The two functions are very close to each other, so this should never happen
    }
    while(overflows_2 > overflows_1){
      scnd_hypergeom_real/=ABS_LIMIT;
      scnd_hypergeom_imag/=ABS_LIMIT;
      //The two functions are very close to each other, so this should never happen
    }
      
		//Calculate prefactor including t^l
		double exp_factor = exp(lnz_fac3*(l+nu_real*0.5)+log(t)*l-overflows_1*log(ABS_LIMIT));
		double exp_phase = lnz_fac3*nu_imag*0.5;
		double exp_real = exp_factor*cos(exp_phase);
		double exp_imag = exp_factor*sin(exp_phase);
		
		double tot_pref_real = MATH_PI_3_2*exp_real;
		double tot_pref_imag = MATH_PI_3_2*exp_imag;
		
		//Calculate gamma factors (large nu -> lnGamma)
		double a_0 = 0.0; double b_0 = 0.0;
		ln_gamma_complex_2(1.5-nu_real*0.5, -nu_imag*0.5, &a_0, &b_0);
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		ln_gamma_complex_2(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);

		exp_factor = exp(a_1+a_2 - a_3 -a_0);
		exp_real = exp_factor*cos(b_1+b_2 - b_3-b_0);
		exp_imag = exp_factor*sin(b_1+b_2 - b_3-b_0);
		
		//Calculate first prefactor
		double pref_first_real = exp_real;
		double pref_first_imag = exp_imag;
		
		//Calculate prefactor of second hypergeom
		exp_factor = pow(z_fac2, nu_real-2.0);
		double temp = log(z_fac2)*nu_imag;
		double pref_scnd_real = exp_factor*cos(temp);
		double pref_scnd_imag = exp_factor*sin(temp);
		
		double c_0 = 0.0; double d_0 = 0.0;
		ln_gamma_complex_2(nu_real*0.5 - 1, nu_imag*0.5, &c_0, &d_0);
		exp_factor = exp(c_0-a_0);
		exp_real = exp_factor*cos(d_0-b_0);
		exp_imag = exp_factor*sin(d_0-b_0);
		
		//Calculate second prefactor
		double pref_num_gamma_real = exp_real;
		double pref_num_gamma_imag = exp_imag;
		double tot_pref_scnd_real = pref_scnd_real*pref_num_gamma_real - pref_scnd_imag*pref_num_gamma_imag;
		double tot_pref_scnd_imag = pref_scnd_real*pref_num_gamma_imag + pref_scnd_imag*pref_num_gamma_real;
		
		//Calculate sum of hypergeoms with their prefactors
		double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*tot_pref_scnd_real - scnd_hypergeom_imag*tot_pref_scnd_imag);
		double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*tot_pref_scnd_imag + scnd_hypergeom_imag*tot_pref_scnd_real);
		
		//Calculate total
		*res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
		*res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
		return;	
	}
}*/
/**
 * The I_l(nu,t) for l = 0 is known analytically. The below implements the corresponding formula
 * */
#define N_BESSEL_L0_COEFFS 16
const double bessel_integral_l0_coeffs[N_BESSEL_L0_COEFFS]={4.,2./3.,1./30., 1.0/1260.,1./90720.,1./9979200.,1./1556755200.,1./326918592000.,1./88921857024000.,1./30411275102208000.,1./12772735542927360000.,1./6463004184721244160000.,1./3877802510832746496000000.,1./2722217362604588040192000000.,1./2210440498434925488635904000000.,1./2055709663544480704431390720000000.};
void bessel_integral_l0(double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
  if( t == 1 ){
    bessel_analytic_limit_atz1(0,nu_real,nu_imag,res_real,res_imag);
    return;
  }
  else if(t<T_MIN_TAYLOR){
    /**
     * 
     *  coeff[i]*(2*pi*cos(nu*pi/2)/2*gamma(nu-2)*(nu-2)/2*gamma(nu-2+2*i)/gamma(nu-2) * (2*i+nu)
     * Even the coefficients (after being multiplied with large gamma functions) decline
     * As such, the series would theoretically converge even towards t=1,
     * but there only very VERY slowly. We expect in that regime the (simple) formula below to hold better anyway
     * 
     * Thus we here explicitly use this formula onyl for t<0.1 
     * 
     * */
    double pref_real,pref_imag; 
    //Use taylor series expansion for small t to avoid numerical cancelations
    // t^20 is already at least of relative accuracy 1e-20, which is less than double
    // However, due to numerical roundoffs, it does not quite reach as low
    if(nu_imag<NU_IMAG_LIMIT){
      double cos_real,cos_imag;
      double gamma_real,gamma_imag;
      cos_real = cos(MATH_PI_HALF*nu_real)*cosh(MATH_PI_HALF*nu_imag);
      cos_imag = -sin(MATH_PI_HALF*nu_real)*sinh(MATH_PI_HALF*nu_imag);
      gamma_complex(nu_real-2.0,nu_imag,&gamma_real,&gamma_imag);
      pref_real = MATH_2_PI*(cos_real*gamma_real-cos_imag*gamma_imag);
      pref_imag = MATH_2_PI*(cos_real*gamma_imag+cos_imag*gamma_real);
    }
    else{
      double ln_cos_real,ln_cos_imag;
      double ln_gam_real,ln_gam_imag;
      double exp_factor;
      ln_cos_complex(MATH_PI_HALF*nu_real,MATH_PI_HALF*nu_imag,&ln_cos_real,&ln_cos_imag);
      ln_gamma_complex_2(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    double tot_pref_real = (nu_real*0.5-1.0)*pref_real-nu_imag*pref_imag*0.5;
    double tot_pref_imag = (nu_real*0.5-1.0)*pref_imag+nu_imag*pref_real*0.5;
    
    double cur_t = 1.0;
    double cur_nu_fac_real,cur_nu_fac_imag;
    double prod_nu_fac_real=1.0;
    double prod_nu_fac_imag=0.0;
    double nu_pref_real;
    double nu_pref_imag;
    double cur_real,cur_imag;
    double sum_real = 0.0;
    double sum_imag = 0.0;
    double temp_real,temp_imag;
    int i;
    for(i=0;i<N_BESSEL_L0_COEFFS;++i){
      nu_pref_real = prod_nu_fac_real;
      nu_pref_imag = prod_nu_fac_imag;
      cur_real = -bessel_integral_l0_coeffs[i]*cur_t*(nu_pref_real*tot_pref_real-nu_pref_imag*tot_pref_imag);
      cur_imag = -bessel_integral_l0_coeffs[i]*cur_t*(nu_pref_real*tot_pref_imag+nu_pref_imag*tot_pref_real);
      
      cur_t *= t*t;
      cur_nu_fac_real = ((nu_real+2.0*i-1.0)*(nu_real+2.0*i)-nu_imag*nu_imag);
      cur_nu_fac_imag = nu_imag*(2.0*nu_real+4.0*i-1.0);
      temp_real = cur_nu_fac_real*prod_nu_fac_real-cur_nu_fac_imag*prod_nu_fac_imag;
      temp_imag = cur_nu_fac_real*prod_nu_fac_imag+cur_nu_fac_imag*prod_nu_fac_real;
      prod_nu_fac_real = temp_real;
      prod_nu_fac_imag = temp_imag;
      //printf("Term %i = %.10e+%.10ej \n",i,cur_real,cur_imag);
      sum_real += cur_real;
      sum_imag += cur_imag;
      //printf("Sum %i = %.10e+%.10ej \n",i,sum_real,sum_imag);
      
    }
    //printf("AT t=%.20e , nu = %.20e+%.20ej , res = %.10e+%.10ej",
    //  t,nu_real,nu_imag,sum_real,sum_imag);
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
  else{
    double log1pt = log(1+t);
    double log1mt = log(1-t);
    double pref_real,pref_imag;
    if(nu_imag<NU_IMAG_LIMIT){
      double cos_real,cos_imag;
      double gamma_real,gamma_imag;
      cos_real = cos(MATH_PI_HALF*nu_real)*cosh(MATH_PI_HALF*nu_imag);
      cos_imag = -sin(MATH_PI_HALF*nu_real)*sinh(MATH_PI_HALF*nu_imag);
      gamma_complex(nu_real-2.0,nu_imag,&gamma_real,&gamma_imag);
      pref_real = MATH_2_PI*(cos_real*gamma_real-cos_imag*gamma_imag);
      pref_imag = MATH_2_PI*(cos_real*gamma_imag+cos_imag*gamma_real);
    }
    else{
      double ln_cos_real,ln_cos_imag;
      double ln_gam_real,ln_gam_imag;
      double exp_factor;
      ln_cos_complex(MATH_PI_HALF*nu_real,MATH_PI_HALF*nu_imag,&ln_cos_real,&ln_cos_imag);
      ln_gamma_complex_2(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    double t_fac_first = exp(log1pt*(2.0-nu_real));
    double t_fac_second = -exp(log1mt*(2.0-nu_real));
    double t_fac_real = t_fac_first*cos(log1pt*nu_imag)+t_fac_second*cos(log1mt*nu_imag);
    double t_fac_imag = -t_fac_first*sin(log1pt*nu_imag)-t_fac_second*sin(log1mt*nu_imag);
    *res_real = (t_fac_real*pref_real-t_fac_imag*pref_imag)/t;
    *res_imag = (t_fac_real*pref_imag+t_fac_imag*pref_real)/t;
    return;
  }
}
/**
 * The I_l(nu,t) for l = 1 is known analytically, the below implements that formula
 * */
#define N_BESSEL_L1_COEFFS 11
const double bessel_integral_l1_coeffs[N_BESSEL_L1_COEFFS]={4./3.,2./15.,1./210., 1.0/11340.,1./997920.,1./129729600.,1./23351328000.,1./5557616064000.,1./1689515283456000.,1./638636777146368000.,1./293772917487329280000.};
//Interestingly, these are related to the l0 coefficients by 1/(2n+1) (with n=index_in_array_starting_from_0+1)
void bessel_integral_l1(double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
  if( t == 1 ){
    bessel_analytic_limit_atz1(1,nu_real,nu_imag,res_real,res_imag);
    return;
  }
  else if(t<T_MIN_TAYLOR){
    /**
     * 
     *  coeff[i]*(2*pi*cos(nu*pi/2)*gamma(nu-2)*(nu-2)/2*gamma(nu-2+2*i)/gamma(nu-2) * (2*i+nu)
     * Even the coefficients (after being multiplied with large gamma functions) decline
     * As such, the series would theoretically converge even towards t=1,
     * but there only very VERY slowly. We expect in that regime the (simple) formula below to hold better anyway
     * 
     * Thus we here explicitly use this formula onyl for t<0.1 
     * 
     * */
    double pref_real,pref_imag; 
    //Use taylor series expansion for small t to avoid numerical cancelations
    // t^21 is already at least of relative accuracy 1e-21, which is less than double
    // However, due to numerical roundoffs, it does not quite reach as low
    if(nu_imag<NU_IMAG_LIMIT){
      double cos_real,cos_imag;
      double gamma_real,gamma_imag;
      cos_real = cos(MATH_PI_HALF*nu_real)*cosh(MATH_PI_HALF*nu_imag);
      cos_imag = -sin(MATH_PI_HALF*nu_real)*sinh(MATH_PI_HALF*nu_imag);
      gamma_complex(nu_real-2.0,nu_imag,&gamma_real,&gamma_imag);
      pref_real = MATH_2_PI*(cos_real*gamma_real-cos_imag*gamma_imag);
      pref_imag = MATH_2_PI*(cos_real*gamma_imag+cos_imag*gamma_real);
    }
    else{
      double ln_cos_real,ln_cos_imag;
      double ln_gam_real,ln_gam_imag;
      double exp_factor;
      ln_cos_complex(MATH_PI_HALF*nu_real,MATH_PI_HALF*nu_imag,&ln_cos_real,&ln_cos_imag);
      ln_gamma_complex_2(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    double tot_pref_real = (nu_real*0.5-1.0)*pref_real-nu_imag*pref_imag*0.5;
    double tot_pref_imag = (nu_real*0.5-1.0)*pref_imag+nu_imag*pref_real*0.5;
    
    double cur_t = t;
    double cur_nu_fac_real,cur_nu_fac_imag;
    double prod_nu_fac_real=1.0;
    double prod_nu_fac_imag=0.0;
    double nu_pref_real;
    double nu_pref_imag;
    double cur_real,cur_imag;
    double sum_real = 0.0;
    double sum_imag = 0.0;
    double temp_real,temp_imag;
    int i;
    for(i=0;i<N_BESSEL_L1_COEFFS;++i){
      nu_pref_real = (nu_real+2.0*i)*prod_nu_fac_real-nu_imag*prod_nu_fac_imag;
      nu_pref_imag = (nu_real+2.0*i)*prod_nu_fac_imag+nu_imag*prod_nu_fac_real;
      cur_real = -bessel_integral_l1_coeffs[i]*cur_t*(nu_pref_real*tot_pref_real-nu_pref_imag*tot_pref_imag);
      cur_imag = -bessel_integral_l1_coeffs[i]*cur_t*(nu_pref_real*tot_pref_imag+nu_pref_imag*tot_pref_real);
      
      cur_t *= t*t;
      cur_nu_fac_real = ((nu_real+2.0*i-1.0)*(nu_real+2.0*i)-nu_imag*nu_imag);
      cur_nu_fac_imag = nu_imag*(2.0*nu_real+4.0*i-1.0);
      temp_real = cur_nu_fac_real*prod_nu_fac_real-cur_nu_fac_imag*prod_nu_fac_imag;
      temp_imag = cur_nu_fac_real*prod_nu_fac_imag+cur_nu_fac_imag*prod_nu_fac_real;
      prod_nu_fac_real = temp_real;
      prod_nu_fac_imag = temp_imag;
      //printf("Term %i = %.10e+%.10ej \n",i,cur_real,cur_imag);
      sum_real += cur_real;
      sum_imag += cur_imag;
      //printf("Sum %i = %.10e+%.10ej \n",i,sum_real,sum_imag);
      
    }
    //printf("AT t=%.20e , nu = %.20e+%.20ej , res = %.10e+%.10ej",
    //  t,nu_real,nu_imag,sum_real,sum_imag);
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
  else{
    double log1pt = log(1+t);
    double log1mt = log(1-t);
    double pref_real,pref_imag;
    if(nu_imag<NU_IMAG_LIMIT){
      double cos_real,cos_imag;
      double gamma_real,gamma_imag;
      cos_real = cos(MATH_PI_HALF*nu_real)*cosh(MATH_PI_HALF*nu_imag);
      cos_imag = -sin(MATH_PI_HALF*nu_real)*sinh(MATH_PI_HALF*nu_imag);
      gamma_complex(nu_real-2.0,nu_imag,&gamma_real,&gamma_imag);
      pref_real = MATH_2_PI*(cos_real*gamma_real-cos_imag*gamma_imag);
      pref_imag = MATH_2_PI*(cos_real*gamma_imag+cos_imag*gamma_real);
      //if(t<1e-3){printf("Pref = %.10e+%.10ej \n",pref_real,pref_imag);}
    }
    else{
      double ln_cos_real,ln_cos_imag;
      double ln_gam_real,ln_gam_imag;
      double exp_factor;
      ln_cos_complex(MATH_PI_HALF*nu_real,MATH_PI_HALF*nu_imag,&ln_cos_real,&ln_cos_imag);
      ln_gamma_complex_2(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    double nu_fac_abs_squared = (4.0-nu_real)*(4.0-nu_real)+nu_imag*nu_imag;
    double tot_pref_real = ((4.0-nu_real)*pref_real-nu_imag*pref_imag)/nu_fac_abs_squared;
    double tot_pref_imag = ((4.0-nu_real)*pref_imag+nu_imag*pref_real)/nu_fac_abs_squared;
    //printf("tot_pref = %.10e+%.10ej \n",tot_pref_real,tot_pref_imag);
    //if(t<1e-3){printf("tot_pref = %.10e+%.10ej \n",tot_pref_real,tot_pref_imag);}
    
    double t_fac_first = exp(log1pt*(2.0-nu_real));
    double t_fac_second = -exp(log1mt*(2.0-nu_real));
    double t_fac_bracket_first_real = (1-t)*(1-t)+nu_real*t;
    double t_fac_bracket_second_real = (1+t)*(1+t)-nu_real*t;
    double first_real = t_fac_first*(cos(log1pt*nu_imag)*t_fac_bracket_first_real+sin(log1pt*nu_imag)*nu_imag*t);
    double first_imag = t_fac_first*(-sin(log1pt*nu_imag)*t_fac_bracket_first_real+cos(log1pt*nu_imag)*nu_imag*t);
    double second_real = t_fac_second*(cos(log1mt*nu_imag)*t_fac_bracket_second_real-sin(log1mt*nu_imag)*nu_imag*t);
    double second_imag = t_fac_second*(-sin(log1mt*nu_imag)*t_fac_bracket_second_real-cos(log1mt*nu_imag)*nu_imag*t);
    
    //if(t<1e-3){printf("t_fac_bracket_first = %.10e \n",t_fac_bracket_first_real);}
    //if(t<1e-3){printf("t_fac_bracket_second = %.10e \n",t_fac_bracket_second_real);}
    //if(t<1e-3){printf("t_fac_first = %.10e \n",t_fac_first);}
    //if(t<1e-3){printf("t_fac_second = %.10e \n",t_fac_second);}
    //if(t<1e-3){printf("first = %.10e+%.10ej \n",first_real,first_imag);}
    //if(t<1e-3){printf("second = %.10e+%.10ej \n",second_real,second_imag);}
    *res_real = ((first_real+second_real)*tot_pref_real-(first_imag+second_imag)*tot_pref_imag)/t/t;
    *res_imag = ((first_real+second_real)*tot_pref_imag+(first_imag+second_imag)*tot_pref_real)/t/t;
    return;
  }
}
#define N_BESSEL_LTAYLOR_COEFFS 40
void bessel_taylor_for_small_t(double l,double nu_real,double nu_imag,double t, double* res_real,double* res_imag){
  double pref_t = pow(t,l);
  double pref_real,pref_imag;
  if(nu_imag>NU_IMAG_LIMIT){
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex_2(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		
		//Divide result by overflow-safe Gamma(l+3/2)
		double exp_factor;
		if(l < BESSEL_INTEGRAL_MAX_L){
			double scnd_den_gamma = gamma_real(l + 1.5);
			exp_factor = exp(lgamma_num_real - a_0+log(t)*l)/scnd_den_gamma;
			
			pref_real = exp_factor*cos(lgamma_num_imag - b_0);
			pref_imag = exp_factor*sin(lgamma_num_imag - b_0);
		}else{
			double scnd_den_gamma_ln = ln_gamma_real(l + 1.5);
			exp_factor = exp(lgamma_num_real - a_0 - scnd_den_gamma_ln+log(t)*l);
			
			pref_real = exp_factor*cos(lgamma_num_imag - b_0);
			pref_imag = exp_factor*sin(lgamma_num_imag - b_0);
		}
	}
	else if (l < BESSEL_INTEGRAL_MAX_L){
		//Calculate Gamma((3-nu)/2) and Gamma(l+3/2) as the denominator
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;
		double scnd_den_gamma = gamma_real(l + 1.5);
		
		//Calculate the numerator
		double gamma_num_real = 0.0; double gamma_num_imag = 0.0;
		gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &gamma_num_real, &gamma_num_imag);
		
		//Calculate fraction of numerator/denominator and multiply with prefactor
		pref_real = (first_den_gamma_real*gamma_num_real+first_den_gamma_imag*gamma_num_imag)/ first_den_gamma_abs_squared / scnd_den_gamma;
		pref_imag = (first_den_gamma_real*gamma_num_imag-first_den_gamma_imag*gamma_num_real)/ first_den_gamma_abs_squared / scnd_den_gamma;	
	}
	else{
		//Only Gamma((3-nu)/2) can be calculated directly
		double first_den_gamma_real = 0.0; double first_den_gamma_imag = 0.0;
		gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &first_den_gamma_real, &first_den_gamma_imag);
		double first_den_gamma_abs_squared = first_den_gamma_real*first_den_gamma_real + first_den_gamma_imag*first_den_gamma_imag;

		//Calculate Gamma(l+nu/2)/Gamma(l+3/2)
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex_2(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		double exp_factor = exp(lgamma_num_real - scnd_den_gamma); 
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);
		
		//Calculate fraction of numerator and denominator
    pref_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		pref_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
	}
  //printf("pref_real = %.10e+%.10ej \n",pref_real,pref_imag);
  
  //Now comes another factor of 2^nu
  double exp_factor = exp(LOG_2*nu_real);
  double tot_pref_real = MATH_PI*MATH_PI*exp_factor*(pref_real*cos(LOG_2*nu_imag)-pref_imag*sin(LOG_2*nu_imag));
  double tot_pref_imag = MATH_PI*MATH_PI*exp_factor*(pref_real*sin(LOG_2*nu_imag)+pref_imag*cos(LOG_2*nu_imag));
  //printf("tot_pref_real = %.10e+%.10ej \n",tot_pref_real,tot_pref_imag);
  //Now comes the actual taylor expansion
  double cur_t = 1.0;
  double cur_fac_real,cur_fac_imag;
  double prod_fac_real = 1.0;
  double prod_fac_imag = 0.0;
  double cur_val_real,cur_val_imag;
  double cur_sum_real=0.0;
  double cur_sum_imag=0.0;
  double temp_real,temp_imag; 
  double taylor_coefficient = 0.5;
  int i;
  for(i=0;i<N_BESSEL_LTAYLOR_COEFFS;++i){
    taylor_coefficient = 1.0/(gamma_real(i+1.0)*pow(2.0,i+1.0));//
    //printf("At i = %i , gamma = %.10e, pow =%.10e ,taylor = %.10e \n",
    //      i,gamma_real(i+1.0),pow(2.0,i+1.0),taylor_coefficient);
    
    cur_val_real = taylor_coefficient*prod_fac_real*cur_t;
    cur_val_imag = taylor_coefficient*prod_fac_imag*cur_t;
    cur_sum_real += cur_val_real;
    cur_sum_imag += cur_val_imag;
    /*printf("cur_val = %.10e + %.10ej :: cur_sum  =%.10e+%.10ej \n",cur_val_real,cur_val_imag,
                                      pref_t*(tot_pref_real*cur_sum_real-tot_pref_imag*cur_sum_imag)
                                      ,pref_t*(tot_pref_real*cur_sum_imag+tot_pref_imag*cur_sum_real)
                                      );*/
    
    cur_t*=t*t;
    //printf("%i prod_fac(old) = %.10e+%.10ej \n",
    //  i,prod_fac_real,prod_fac_imag);
      
    cur_fac_real = ((-1.0+nu_real+2.0*i)*(l+0.5*nu_real+i)-0.5*nu_imag*nu_imag)/(l+1.5+i);
    cur_fac_imag = ((-1.0+nu_real+2.0*i)*0.5*nu_imag+(l+0.5*nu_real+i)*nu_imag)/(l+1.5+i);
    
    temp_real = (cur_fac_real*prod_fac_real-cur_fac_imag*prod_fac_imag);
    temp_imag = (cur_fac_imag*prod_fac_real+cur_fac_real*prod_fac_imag);
    prod_fac_real = temp_real;
    prod_fac_imag = temp_imag;
    
  }
  //printf("t= %.10e, pref_t = %.10e, tot_pref = %.10e+%.10ej,cur_sum = %.10e+%.10ej \n",
  //      t,pref_t,tot_pref_real,tot_pref_imag,cur_sum_real,cur_sum_imag);
  *res_real = pref_t*(tot_pref_real*cur_sum_real-tot_pref_imag*cur_sum_imag);
  *res_imag = pref_t*(tot_pref_imag*cur_sum_real+tot_pref_real*cur_sum_imag);
  return;
}
/**
 * This function calculates the integral over the bessel functions.
 * 
 * When getting into regimes of t<<1 or (1-t*t)<<1, we have to switch the method of calculation
 * This is done using a BESSEL_INTEGRAL_T_SWITCH which decides how high t² should be for a switch to occur
 * (The bessel_integral_transform version should be preferred)
 * */
void bessel_integral(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	if(t==1){
		bessel_analytic_limit_atz1(l,nu_real,nu_imag,res_real,res_imag);
		return;
	}
	else if (t*t < BESSEL_INTEGRAL_T_SWITCH ||(nu_real==0.0 && nu_imag==0.0)){
		bessel_integral_lowt_lownu(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
	else{
		bessel_integral_hight_lownu(l, nu_real, nu_imag, t, res_real, res_imag);
		return;
	}
}
/**
 * This function calculates the integral over the bessel functions using the transformed versions of the equations
 * 
 * When l=0 or l=1 the function is a very simple analytical one
 * 
 * When getting into regimes of t<<1 or (1-t*t)<<1, we have to switch the method of calculation
 * This is done using a BESSEL_INTEGRAL_T_SWITCH which decides how high t² should be for a switch to occur
 * */
void bessel_integral_transform(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	if( l == 0){
    bessel_integral_l0(nu_real,nu_imag,t,res_real,res_imag);
    return;
  }
  if( l == 1){
    bessel_integral_l1(nu_real,nu_imag,t,res_real,res_imag);
    return;
  }
  if (t == 1){
		bessel_analytic_limit_atz1(l,nu_real,nu_imag,res_real,res_imag);
		return;
	}
  if (t<T_MIN_TAYLOR){
    bessel_taylor_for_small_t(l,nu_real,nu_imag,t,res_real,res_imag);
  }
	else if(t*t<BESSEL_INTEGRAL_T_SWITCH_SAFE){
    //printf("low t = %.10e \n",t);
		//bessel_integral_lowt_lownu(l, nu_real, nu_imag, t, res_real,res_imag); // Less efficient, but sometimes useful (?)
		bessel_integral_lowt_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
  else if(t*t<BESSEL_INTEGRAL_T_SWITCH){
    //printf("Safe t = %.10e \n",t);
    bessel_integral_lowt_transform_safe(l, nu_real, nu_imag, t, res_real,res_imag);
    return;
  }
	else{
    //printf("high t = %.10e \n",t);
		bessel_integral_hight_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		//bessel_integral_hight_transform_safe(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
}
/**
 * This function calculates the integral over the bessel functions using the transformed versions of the equations
 * 
 * When l=0 or l=1 the function is a very simple analytical one
 * 
 * When getting into regimes of t<<1 or (1-t*t)<<1, we have to switch the method of calculation
 * This is done using a BESSEL_INTEGRAL_T_SWITCH which decides how high t² should be for a switch to occur
 * */
void bessel_integral_transform_safety(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag, short* overflow_flag){
	if( l == 0){
    bessel_integral_l0(nu_real,nu_imag,t,res_real,res_imag);
    return;
  }
  if( l == 1){
    bessel_integral_l1(nu_real,nu_imag,t,res_real,res_imag);
    return;
  }
  if (t == 1){
		bessel_analytic_limit_atz1(l,nu_real,nu_imag,res_real,res_imag);
		return;
	}
  if (t<T_MIN_TAYLOR){
    bessel_taylor_for_small_t(l,nu_real,nu_imag,t,res_real,res_imag);
  }
	else if(t*t<BESSEL_INTEGRAL_T_SWITCH_SAFE && !(*overflow_flag)){
    //printf("low t = %.10e \n",t);
		//bessel_integral_lowt_lownu(l, nu_real, nu_imag, t, res_real,res_imag); // Less efficient, but sometimes useful (?)
		bessel_integral_lowt_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
  else if(t*t<BESSEL_INTEGRAL_T_SWITCH){
    //printf("Safe t = %.10e \n",t);
    bessel_integral_lowt_transform_safe_flag(l, nu_real, nu_imag, t, res_real,res_imag,overflow_flag);
    return;
  }
	else{
    //printf("high t = %.10e \n",t);
		bessel_integral_hight_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		//bessel_integral_hight_transform_safe(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
}

#define NUM_HYPERGEOM_TAYLOR_COEFF 15
#define BESSEL_HYPERGEOM_TAYLOR 0.01
void bessel_hypergeom_l0(double nu_real,double nu_imag,double t,double* res_real,double* res_imag){
  if(t>BESSEL_HYPERGEOM_TAYLOR){
    double log1pt = log(1.0+t);
    double log1mt = log(1.0-t); //Can be -inf for t=1.0
    double exp_factor = exp(log1pt*(2-nu_real));
    double first_real = -exp_factor*cos(-log1pt*nu_imag);
    double first_imag = -exp_factor*sin(-log1pt*nu_imag);
    exp_factor = exp(log1mt*(2-nu_real));
    double second_real = exp_factor*cos(-log1mt*nu_imag);
    double second_imag = exp_factor*sin(-log1mt*nu_imag);
    double den_abs_sq = (nu_real-2.0)*(nu_real-2.0)+nu_imag*nu_imag;
    
    *res_real = 0.5/t*((first_real+second_real)*(nu_real-2.0)+(first_imag+second_imag)*nu_imag)/den_abs_sq;
    *res_imag = 0.5/t*((first_imag+second_imag)*(nu_real-2.0)-(first_real+second_real)*nu_imag)/den_abs_sq;
    return;
  }
  else{
    //If t is too small the above depends on precise cancellations
    //Instead, we use this below, which expands in small t
    double cur_nu_fac_real=1.0,cur_nu_fac_imag=0.0;
    double cur_t_fac = 1.0;
    double sum_real = 0.0;
    double sum_imag = 0.0;
    double temp_real,temp_imag;
    double temp_nu_real,temp_nu_imag;
    int i;
    for(i=0;i<NUM_HYPERGEOM_TAYLOR_COEFF;++i){
      sum_real += cur_nu_fac_real*cur_t_fac/gamma_real(2*i+2);
      sum_imag += cur_nu_fac_imag*cur_t_fac/gamma_real(2*i+2);
      
      cur_t_fac*=t*t;
      temp_nu_real = (nu_real+2*i-1)*(nu_real+2*i)-nu_imag*nu_imag;
      temp_nu_imag = (nu_real+2*i-1)*nu_imag+(nu_real+2*i)*nu_imag;
      temp_real = (temp_nu_real*cur_nu_fac_real-temp_nu_imag*cur_nu_fac_imag);
      temp_imag = (temp_nu_real*cur_nu_fac_imag+temp_nu_imag*cur_nu_fac_real);
      cur_nu_fac_real = temp_real;
      cur_nu_fac_imag = temp_imag;
    }
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
}
void bessel_hypergeom_l0_p1(double nu_real,double nu_imag,double t,double* res_real,double* res_imag){
  if(t>BESSEL_HYPERGEOM_TAYLOR){
    double log1pt = log(1.0+t);
    double log1mt = log(1.0-t); //Can be -inf for t=1.0
    double den_abs_sq = (nu_real-2.0)*(nu_real-2.0)+nu_imag*nu_imag;
    
    double exp_factor = exp(log1mt*(1-nu_real));
    double first_real = exp_factor*cos(-log1mt*nu_imag);
    double first_imag = exp_factor*sin(-log1mt*nu_imag);
    exp_factor = exp(log1pt*(1-nu_real));
    double second_real = exp_factor*cos(-log1pt*nu_imag);
    double second_imag = exp_factor*sin(-log1pt*nu_imag);
    
    double first_bracket_real = (1-t)*first_real-(1+t)*second_real;
    double first_bracket_imag = (1-t)*first_imag-(1+t)*second_imag;
    double pref_first_bracket_real = 1/t*((nu_real-1.0)*(nu_real-2.0)+nu_imag*nu_imag)/den_abs_sq;
    double pref_first_bracket_imag = 1/t*(nu_imag*(nu_real-2.0)-(nu_real-1.0)*nu_imag)/den_abs_sq;
    
    double tot_first_bracket_real = first_bracket_real*pref_first_bracket_real-first_bracket_imag*pref_first_bracket_imag;
    double tot_first_bracket_imag = first_bracket_imag*pref_first_bracket_real+first_bracket_real*pref_first_bracket_imag;
    
    double second_bracket_real = first_real+second_real;
    double second_bracket_imag = first_imag+second_imag;
    
    double tot_brackets_real = tot_first_bracket_real+second_bracket_real;
    double tot_brackets_imag = tot_first_bracket_imag+second_bracket_imag;
    
    den_abs_sq = nu_real*nu_real+nu_imag*nu_imag;
    *res_real = 0.5*(tot_brackets_real*nu_real+tot_brackets_imag*nu_imag)/den_abs_sq;
    *res_imag = 0.5*(tot_brackets_imag*nu_real-tot_brackets_real*nu_imag)/den_abs_sq;
    return;
  }
  else{
    //If t is too small the above depends on precise cancellations
    //Instead, we use this below, which expands in small t
    double cur_nu_fac_real=1.0,cur_nu_fac_imag=0.0;
    double cur_t_fac = 1.0;
    double sum_real = 0.0;
    double sum_imag = 0.0;
    double temp_real,temp_imag;
    double temp_nu_real,temp_nu_imag;
    int i;
    for(i=0;i<NUM_HYPERGEOM_TAYLOR_COEFF;++i){
      sum_real += cur_nu_fac_real*cur_t_fac/gamma_real(2*i+2);
      sum_imag += cur_nu_fac_imag*cur_t_fac/gamma_real(2*i+2);
      
      cur_t_fac*=t*t;
      temp_nu_real = (nu_real+2*i-1)*(nu_real+2*i+2)-nu_imag*nu_imag;
      temp_nu_imag = (nu_real+2*i-1)*nu_imag+(nu_real+2*i+2)*nu_imag;
      temp_real = (temp_nu_real*cur_nu_fac_real-temp_nu_imag*cur_nu_fac_imag);
      temp_imag = (temp_nu_real*cur_nu_fac_imag+temp_nu_imag*cur_nu_fac_real);
      cur_nu_fac_real = temp_real;
      cur_nu_fac_imag = temp_imag;
    }
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
}
#define GIANT_VAL 1e120
#define NON_INVERTIBLE_LIMIT 1e-25
//1e-20
int bessel_integral_recursion_complicated(int l_max,int l_recursion_max,double nu_real, double nu_imag, double t,double bi_allowed_error,double* bi_real,double* bi_imag,double* max_t,double* initial_abs,ErrorMsg errmsg){
  
  //clock_t start = clock();
  double z = t*t;
  double a_real = 0.5*(nu_real-1.0);
  double a_imag = 0.5*nu_imag;
  double m11_real = 1.0,m11_imag=0.0;
  double m12_real = 0.0,m12_imag=0.0;
  double m21_real = 0.0,m21_imag=0.0;
  double m22_real = 1.0,m22_imag=0.0;
  double f000_real=1.0,f000_imag=0.0;
  double f010_real=1.0,f010_imag=0.0;
  double temp1_real,temp1_imag;
  double temp2_real,temp2_imag;
  double temp3_real,temp3_imag;
  double temp4_real,temp4_imag;
  double tempf_real,tempf_imag;
  bi_real[l_max] = 1.0;
  bi_imag[l_max] = 0.0;
  int index_l;
  for(index_l=l_recursion_max;index_l>0;--index_l){
    double c = 1.5+index_l-1;
    double alpha_real = a_real/c;
    double alpha_imag = a_imag/c;
    double z_alpha_real = z-alpha_real;
    double o_alpha_real = 1.+alpha_real;
    //Calculate new matrix
    //printf("a = %.10e+%.10ej , c = %.10e , z = %.10e \n",a_real,a_imag,c,z);
    
    //new 11
    temp1_real = (z_alpha_real)*m11_real+alpha_imag*m11_imag+(o_alpha_real)*(1-z)*m21_real-alpha_imag*(1-z)*m21_imag;
    temp1_imag = (z_alpha_real)*m11_imag-alpha_imag*m11_real+(o_alpha_real)*(1-z)*m21_imag+alpha_imag*(1-z)*m21_real;
    
    //new 12
    temp2_real = (z_alpha_real)*m12_real+alpha_imag*m12_imag+(o_alpha_real)*(1-z)*m22_real-alpha_imag*(1-z)*m22_imag;
    temp2_imag = (z_alpha_real)*m12_imag-alpha_imag*m12_real+(o_alpha_real)*(1-z)*m22_imag+alpha_imag*(1-z)*m22_real;
    //new 21
    temp3_real = -alpha_real*m11_real+alpha_imag*m11_imag+(o_alpha_real)*m21_real-alpha_imag*m21_imag;
    temp3_imag = -alpha_real*m11_imag-alpha_imag*m11_real+(o_alpha_real)*m21_imag+alpha_imag*m21_real;
    
    //new 22
    temp4_real = -alpha_real*m12_real+alpha_imag*m12_imag+(o_alpha_real)*m22_real-alpha_imag*m22_imag;
    temp4_imag = -alpha_real*m12_imag-alpha_imag*m12_real+(o_alpha_real)*m22_imag+alpha_imag*m22_real;
    
    //Advance (f000,f010) vector
    tempf_real = (z_alpha_real)*f000_real+alpha_imag*f000_imag+(o_alpha_real)*(1-z)*f010_real-alpha_imag*(1-z)*f010_imag;
    tempf_imag = (z_alpha_real)*f000_imag-alpha_imag*f000_real+(o_alpha_real)*(1-z)*f010_imag+alpha_imag*(1-z)*f010_real;
    //Save f000 results and obtain f010 results
    bi_real[index_l-1] = tempf_real;
    bi_imag[index_l-1] = tempf_imag;
    //printf("f = %+.10e + %+.10ej  -- %+.10e + %+.10ej \n",f000_real,f000_imag,f010_real,f010_imag);
    tempf_real = -alpha_real*f000_real+alpha_imag*f000_imag+(o_alpha_real)*f010_real-alpha_imag*f010_imag;
    tempf_imag = -alpha_real*f000_imag-alpha_imag*f000_real+(o_alpha_real)*f010_imag+alpha_imag*f010_real;
    //Prepare (f000,f010) vector for next step
    f010_real = tempf_real;
    f010_imag = tempf_imag;
    f000_real = bi_real[index_l-1];
    f000_imag = bi_imag[index_l-1];
    //printf("f' = %+.10e + %+.10ej  -- %+.10e + %+.10ej \n",f000_real,f000_imag,f010_real,f010_imag);
    
    //Prepare matrix for next step
    m11_real = temp1_real;
    m11_imag = temp1_imag;
    m12_real = temp2_real;
    m12_imag = temp2_imag;
    m21_real = temp3_real;
    m21_imag = temp3_imag;
    m22_real = temp4_real;
    m22_imag = temp4_imag;
    //printf("%+.10e + %+.10ej  -- %+.10e + %+.10ej \n",m11_real,m11_imag,m12_real,m12_imag);
    //printf("%+.10e + %+.10ej  -- %+.10e + %+.10ej \n",m21_real,m21_imag,m22_real,m22_imag);
    //printf("%+.10e + %+.10ej  -- %+.10e + %+.10ej \n",f000_real,f000_imag,f010_real,f010_imag);
    
  }
  //printf("Found final f000 = %.10e+%.10ej and f010 = %.10e+%.10ej \n",f000_real,f000_imag,f010_real,f010_imag);
  
  //Now we know the backward recursion-vector from (1,1), including all bessel values inbetween
  //and also the transition matrix.
  double ad_real = (m11_real*m22_real-m11_imag*m22_imag);
  double ad_imag = (m11_real*m22_imag+m11_imag*m22_real);
  double bc_real = (m12_real*m21_real-m12_imag*m21_imag);
  double bc_imag = (m12_real*m21_imag+m12_imag*m21_real);
  
  
  if(fabs(ad_real)>GIANT_VAL || fabs(ad_imag)>GIANT_VAL || fabs(bc_real)>GIANT_VAL || fabs(bc_imag)>GIANT_VAL){
    ad_real/=GIANT_VAL;
    ad_imag/=GIANT_VAL;
    bc_real/=GIANT_VAL;
    bc_imag/=GIANT_VAL;
    //This leaves the relative factor invariant
  }
  if(fabs(ad_real)<1./GIANT_VAL || fabs(ad_imag)<1./GIANT_VAL || fabs(bc_real)<1./GIANT_VAL || fabs(bc_imag)<1./GIANT_VAL){
    ad_real*=GIANT_VAL;
    ad_imag*=GIANT_VAL;
    bc_real*=GIANT_VAL;
    bc_imag*=GIANT_VAL;
    //This leaves the relative factor invariant
  }
  double det_abs = (ad_real-bc_real)*(ad_real-bc_real)+(ad_imag-bc_imag)*(ad_imag-bc_imag);
  double ad_abs = ad_real*ad_real+ad_imag*ad_imag;
  double bc_abs = bc_real*bc_real+bc_imag*bc_imag;
  //If this matrix is invertible, we are in the forward recursion case
  //Otherwise, we have already found the results by multiplying with lambda
  double res_real,res_imag;
  bessel_hypergeom_l0(nu_real,nu_imag,t,&res_real,&res_imag);
  //printf("ad = %.10e+%.10ej , bc = %.10e+%.10ej \n",ad_real,ad_imag,bc_real,bc_imag);
  //printf("det = %.10e+%.10ej \n",ad_real-bc_real,ad_imag-bc_imag);
  //printf("Det_abs = %.10e , ad_abs = %.10e , bc_abs = %.10e \n",det_abs,ad_abs,bc_abs);
  ///printf("Limit val = %.10e \n",det_abs/(ad_abs+bc_abs));
  if(det_abs/(ad_abs+bc_abs) < NON_INVERTIBLE_LIMIT){
    
    ///printf("Backward recursion assumed \n");
    double lambda_real = (f000_real*res_real+f000_imag*res_imag)/(res_real*res_real+res_imag*res_imag);
    double lambda_imag = (f000_imag*res_real-f000_real*res_imag)/(res_real*res_real+res_imag*res_imag);
    
    double lambda_abs_sq = lambda_real*lambda_real+lambda_imag*lambda_imag;
    
    //printf("Found lambda = %.10e+%.10ej \n",lambda_real,lambda_imag);
    double temp_real,temp_imag;
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      //printf("Before Bi = %.10e+%.10ej \n",bi_real[index_l],bi_imag[index_l]);
      temp_real = (bi_real[index_l]*lambda_real+bi_imag[index_l]*lambda_imag)/lambda_abs_sq;
      temp_imag = (bi_imag[index_l]*lambda_real-bi_real[index_l]*lambda_imag)/lambda_abs_sq;
      bi_real[index_l] = temp_real;
      bi_imag[index_l] = temp_imag;
      //printf("After Bi = %.10e+%.10ej \n",temp_real,temp_imag);
    }
    double err = sqrt(((bi_real[0]-res_real)*(bi_real[0]-res_real)+(bi_imag[0]-res_imag)*(bi_imag[0]-res_imag))/(res_real*res_real+res_imag*res_imag));
    if(err>bi_allowed_error){
      sprintf(errmsg,"%s(L:%d) : Backwards recursion for bessel integrals returned unnaturally high error. Outside of acceptable parameter region assumed.",__func__,__LINE__);
      return _FAILURE_;
    }
  }
  else{
    ///printf("Forward recursion assumed \n");
    bi_real[0] = res_real;
    bi_imag[0] = res_imag; 
    double temp1_real,temp1_imag;
    double temp2_real,temp2_imag;
    double begin_real,begin_imag;
    bessel_hypergeom_l0_p1(nu_real,nu_imag,t,&begin_real,&begin_imag);
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      double c = 1.5+index_l;
      double alpha_real = a_real/c;
      double alpha_imag = a_imag/c;
      //printf("Alpha = %.10e+%.10ej \n",alpha_real,alpha_imag);
      
      double alpha_m1_abs_sq = (1-alpha_real)*(1-alpha_real)+alpha_imag*alpha_imag;
      double alpha_m1_sq_real = 1.0-alpha_real*alpha_real+alpha_imag*alpha_imag;
      double alpha_m1_sq_imag =    -alpha_real*alpha_imag-alpha_imag*alpha_real;
      //printf("1-alpha^2 = %.10e+%.10ej \n",alpha_m1_sq_real,alpha_m1_sq_imag);
      
      //printf("res = %.10e+%.10ej , begin = %.10e+%.10ej \n",res_real,res_imag,begin_real,begin_imag);
      double den_abs_sq = alpha_m1_sq_real*alpha_m1_sq_real+alpha_m1_sq_imag*alpha_m1_sq_imag;
      double c_real = (alpha_real*alpha_m1_sq_real+alpha_imag*alpha_m1_sq_imag)/den_abs_sq/z;
      double c_imag = (alpha_imag*alpha_m1_sq_real-alpha_real*alpha_m1_sq_imag)/den_abs_sq/z;
      double d_real = ( (z-alpha_real)*alpha_m1_sq_real-alpha_imag*alpha_m1_sq_imag)/den_abs_sq/z;
      double d_imag = (-(z-alpha_real)*alpha_m1_sq_imag-alpha_imag*alpha_m1_sq_real)/den_abs_sq/z;
      //printf("-> %.10e+%.10e , %.10e+%.10e \n",(1-alpha_real)/alpha_m1_abs_sq/z,alpha_imag/alpha_m1_abs_sq/z,-(1-z)*(1-alpha_real)/alpha_m1_abs_sq/z,-(1-z)*alpha_imag/alpha_m1_abs_sq/z);
      //printf("-> %.10e+%.10e , %.10e+%.10e \n \n",c_real,c_imag,d_real,d_imag);
      temp1_real = ((1-alpha_real)*res_real-alpha_imag*res_imag-(1-z)*(1-alpha_real)*begin_real+(1-z)*alpha_imag*begin_imag)/alpha_m1_abs_sq/z;
      temp1_imag = ((1-alpha_real)*res_imag+alpha_imag*res_real-(1-z)*(1-alpha_real)*begin_imag-(1-z)*alpha_imag*begin_real)/alpha_m1_abs_sq/z;
      temp2_real = (c_real*res_real-c_imag*res_imag+d_real*begin_real-d_imag*begin_imag);
      temp2_imag = (c_real*res_imag+c_imag*res_real+d_real*begin_imag+d_imag*begin_real);
      bi_real[index_l+1] = temp1_real;
      bi_imag[index_l+1] = temp1_imag;
      begin_real = temp2_real;
      begin_imag = temp2_imag;
      res_real = bi_real[index_l+1];
      res_imag = bi_imag[index_l+1];
      
      //printf("res = %.10e+%.10ej , begin = %.10e+%.10ej \n",res_real,res_imag,begin_real,begin_imag);
    }
  }
  //The bi_real and bi_imag arrays are now filled with Hypergeometric2F1((nu-1)/2,l+nu/2,l+3/2,t*t),
  // but now we want to include the prefactors
  double temp_real,temp_imag;
  double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1.0));
  double pref_pref_real = exp_factor*cos(LOG_2*nu_imag);
  double pref_pref_imag = exp_factor*sin(LOG_2*nu_imag);
  double start_frac_real,start_frac_imag;
  if(nu_imag<NU_IMAG_LIMIT){
    double den_gamma_real,den_gamma_imag;
    double starting_gamma_real,starting_gamma_imag;
    gamma_complex(1.5-0.5*nu_real,-0.5*nu_imag,&den_gamma_real,&den_gamma_imag);
    gamma_complex(nu_real*0.5,nu_imag*0.5,&starting_gamma_real,&starting_gamma_imag);
    double den_gamma_abs_sq = den_gamma_real*den_gamma_real+den_gamma_imag*den_gamma_imag;
    start_frac_real = (den_gamma_real*starting_gamma_real+den_gamma_imag*starting_gamma_imag)/den_gamma_abs_sq;
    start_frac_imag = (den_gamma_real*starting_gamma_imag-den_gamma_imag*starting_gamma_real)/den_gamma_abs_sq;
  }
  else{
    double a0,b0;
    double a1,b1;
    ln_gamma_complex_2(1.5-0.5*nu_real,-0.5*nu_imag,&a0,&b0);
    ln_gamma_complex_2(0.5*nu_real,0.5*nu_imag,&a1,&b1);
    exp_factor = exp(a1-a0);
    start_frac_real = exp_factor*cos(b1-b0);
    start_frac_imag = exp_factor*sin(b1-b0);
  }
  temp_real = 2.0/sqrt(MATH_PI)*(pref_pref_real*start_frac_real-pref_pref_imag*start_frac_imag);
  temp_imag = 2.0/sqrt(MATH_PI)*(pref_pref_imag*start_frac_real+pref_pref_real*start_frac_imag);
  start_frac_real= temp_real;
  start_frac_imag= temp_imag;
  
  double cur_frac_real = start_frac_real;
  double cur_frac_imag = start_frac_imag;
  for(index_l=0;index_l<=l_max;++index_l){
    temp_real = cur_frac_real*bi_real[index_l]-cur_frac_imag*bi_imag[index_l];
    temp_imag = cur_frac_imag*bi_real[index_l]+cur_frac_real*bi_imag[index_l];
    bi_real[index_l] = temp_real;
    bi_imag[index_l] = temp_imag;
    
    if(bi_real[index_l]*bi_real[index_l]+bi_imag[index_l]*bi_imag[index_l]<BESSEL_EPSILON*initial_abs[index_l]){
      max_t[index_l] = t;
      //printf("Exiting :: %4d \n",index_l);
    }
    temp_real = t*(cur_frac_real*(index_l+0.5*nu_real)-cur_frac_imag*0.5*nu_imag)/(index_l+1.5);
    temp_imag = t*(cur_frac_imag*(index_l+0.5*nu_real)+cur_frac_real*0.5*nu_imag)/(index_l+1.5); 
    cur_frac_real = temp_real;
    cur_frac_imag = temp_imag;
  }
  //clock_t end = clock();
  //printf(" -> Back complicated took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
  return _SUCCESS_;
}
void bessel_integral_recursion_initial_abs(int l_max,double nu_real,double nu_imag,double* abi_real,double* abi_imag,double* initial_abs){
  
  //clock_t start = clock();
  double bi_real,bi_imag,bi_next_real,bi_next_imag;
  bessel_integral_l0(
                     nu_real,
                     nu_imag,
                     1.0,
                     &bi_real,
                     &bi_imag
                     );
  bessel_integral_l1(
                     nu_real,
                     nu_imag,
                     1.0,
                     &bi_next_real,
                     &bi_next_imag
                     );
  double bi_next_next_real;
  double bi_next_next_imag;
  double temp_real, temp_imag;
  int index_l;
  for(index_l=0;index_l<=l_max;++index_l){
    double l = (double)index_l;
    abi_real[index_l]= bi_real;
    abi_imag[index_l]= bi_imag;
    initial_abs[index_l]=bi_real*bi_real+bi_imag*bi_imag;
    //Explicitly t=1 is important for criterium of 'vanishing' function 
    double den_real = (3.0+l-nu_real*0.5);
    double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
    double bi_factor_real = -((l+nu_real*0.5)*bi_real-nu_imag*bi_imag*0.5);
    double bi_factor_imag = -((l+nu_real*0.5)*bi_imag+nu_imag*bi_real*0.5);
    double bi_factor_next_real = 2*(l+1.5)*bi_next_real;
    double bi_factor_next_imag = 2*(l+1.5)*bi_next_imag;
    temp_real = (bi_factor_next_real+bi_factor_real);
    temp_imag = (bi_factor_next_imag+bi_factor_imag);
    bi_next_next_real = (temp_real*den_real-temp_imag*nu_imag*0.5)/den_abs_squared;
    bi_next_next_imag = (temp_imag*den_real+temp_real*nu_imag*0.5)/den_abs_squared;
    bi_real = bi_next_real;
    bi_imag = bi_next_imag;
    bi_next_real = bi_next_next_real;
    bi_next_imag = bi_next_next_imag;
  }
  //clock_t end = clock();
  //printf(" -> Initial Abs took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
  return;
} 
void bessel_integral_recursion_taylor(int l_max,double nu_real,double nu_imag,double t,double* max_t,double* initial_abs,double* bi_real,double* bi_imag){
  //clock_t start = clock();
  double res_real,res_imag;
  int index_l;
  for(index_l=0;index_l<=l_max;++index_l){
    if(t<max_t[index_l]){
      //TODO :: fix ?
      bi_real[index_l]=0.0;
      bi_imag[index_l]=0.0;
    }
    else{
      double l = (double)index_l;
      bessel_integral_transform(l,nu_real,nu_imag,t,&res_real,&res_imag);
      if(res_real*res_real+res_imag*res_imag<BESSEL_EPSILON*initial_abs[index_l]){
        max_t[index_l]=t;
      }
      bi_real[index_l] = res_real;
      bi_imag[index_l] = res_imag;
      //End exit check
    }
    //End max t
  }
  //End l
  //clock_t end = clock();
  //printf(" -> Recursion Taylor took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
}


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
                                            ){
                                              
  //clock_t start = clock();
  double bi_real,bi_imag,bi_next_real,bi_next_imag,bi_next_next_real,bi_next_next_imag;
  double temp_real,temp_imag;
  int bessel_recursion_backward_min_l_step = (int)(nu_fraction*(bessel_recursion_backward_min_l_step_high_nu-bessel_recursion_backward_min_l_step_low_nu)+bessel_recursion_backward_min_l_step_low_nu);
  //if(index_coeff>105 && t>0.945){printf("step = %d",bessel_recursion_backward_min_l_step);}
  double min_factor = bessel_recursion_backward_min_l_step*T_SWITCH_FORWARD_RECURSION_SIMPLE;
  int bessel_recursion_l_step = (int)MIN(min_factor/t,bessel_recursion_backward_max_l_step);
  
  int num_recursion_steps = (int)ceil((double)(l_recursion_max+1)/(double)bessel_recursion_l_step);
  int index_recursion_step;
  for(index_recursion_step = 0 ;index_recursion_step < num_recursion_steps;++index_recursion_step){
    int l_start = MAX(l_recursion_max-index_recursion_step*bessel_recursion_l_step,2);
    //printf("Starting at %4d \n",l_start);
    //while(t<max_t[l_start]){
    //  l_start--;
    //}
    //printf("Current rec_step  = %4d , l_max = %4d , l_start = %4d , dl = %4d \n",index_recursion_step,l_max,l_start,(index_recursion_step+1)*bessel_recursion_l_step);
    bessel_integral_transform_safety(
                       l_start,
                       nu_real,
                       nu_imag,
                       t,
                       &bi_real,
                       &bi_imag,
                       overflow_flag);
    bessel_integral_transform_safety(
                       l_start-1,
                       nu_real,
                       nu_imag,
                       t,
                       &bi_next_real,
                       &bi_next_imag,
                       overflow_flag);
    
    abi_real[l_start] = bi_real;
    abi_imag[l_start] = bi_imag;
    abi_real[l_start-1] = bi_next_real;
    abi_imag[l_start-1] = bi_next_imag;
    int index_l;
    for(index_l=l_start;index_l>=(l_recursion_max-(index_recursion_step+1)*bessel_recursion_l_step) && index_l>=2;--index_l){  

      //printf("Allowed, since max_t[%i]=%.10e \n",index_l,max_t[index_l]);
      double l = (double)index_l;

      //printf("Starting with %.10e+%.10ej (l %2d) and %.10e+%.10ej (l-1 %2d) \n",
      //       bi_real,bi_imag,(int)l_start,bi_next_real,bi_next_imag,(int)l_start-1);
      //Explicitly t=1 is important for criterium of 'vanishing' function 
      
      double den_real = (l-2.0+nu_real*0.5);
      double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
      double bi_factor_real = -((1.0+l-nu_real*0.5)*bi_real+nu_imag*bi_imag*0.5);
      double bi_factor_imag = -((1.0+l-nu_real*0.5)*bi_imag-nu_imag*bi_real*0.5);
      double bi_factor_next_real = (1+t*t)/t*(l-0.5)*bi_next_real;
      double bi_factor_next_imag = (1+t*t)/t*(l-0.5)*bi_next_imag;
      //printf("%4d bi facs = %.10e+%.10ej and %.10e+%.10ej \n",index_l,bi_factor_real,bi_factor_imag,bi_factor_next_real,bi_factor_next_imag);
      
      temp_real = (bi_factor_next_real+bi_factor_real);
      temp_imag = (bi_factor_next_imag+bi_factor_imag);
      bi_next_next_real = (temp_real*den_real+temp_imag*nu_imag*0.5)/den_abs_squared;
      bi_next_next_imag = (temp_imag*den_real-temp_real*nu_imag*0.5)/den_abs_squared;
      bi_real = bi_next_real;
      bi_imag = bi_next_imag;
      bi_next_real = bi_next_next_real;
      bi_next_imag = bi_next_next_imag;
      
      abi_real[index_l-2] = bi_next_next_real;
      abi_imag[index_l-2] = bi_next_next_imag;
      /*if(bessel_recursion_error_check>1){
        bessel_integral_transform(l-2.0,nu_real,nu_imag,t,&res_real,&res_imag);
        double err = sqrt((((bi_next_real-res_real)*(bi_next_real-res_real))
              +((bi_next_imag-res_imag)*(bi_next_imag-res_imag)))
              /(res_real*res_real+res_imag*res_imag));
        
        //printf("Err %.10e :: %.10e+%.10ej vs %.10e+%.10ej \n",
        //       err,bi_next_real,bi_next_imag,res_real,res_imag);
        
        class_test(err>bi_allowed_error,
              errmsg,
              "recursion MIDDLE on backward recursion exceeds maximum allowed error at l = %4d, t = %.10e (err = %.10e, seed offset %1d , last seed = %2d) \n Compare wrong %.10e+%.10ei to wanted %.10e+%.10ei",
              (int)l-2,t,err,index_recursion_step,(int)l_start+2,bi_next_real,bi_next_imag,res_real,res_imag);
      }*/
      if(bi_next_real*bi_next_real+bi_next_imag*bi_next_imag<BESSEL_EPSILON*initial_abs[index_l]){ 
        /*max_t[index_l]=t;
        if(index_l==*l_max_cur){
          (*l_max_cur)--;
        }*/
        /*if(bessel_recursion_error_check>0){
          bessel_integral_transform(l-2.0,nu_real,nu_imag,t,&res_real,&res_imag);
          double err = sqrt((((bi_next_real-res_real)*(bi_next_real-res_real))
              +((bi_next_imag-res_imag)*(bi_next_imag-res_imag)))
              /(res_real*res_real+res_imag*res_imag));
          class_test(err>bi_allowed_error,
              errmsg,
              "recursion exit on backward recursion exceeds maximum allowed error at l = %4d, t = %.10e (err = %.10e, seed offset %1d , last seed = %2d) \n Compare wrong %.10e+%.10ei to wanted %.10e+%.10ei",
              (int)l-2,t,err,index_recursion_step,(int)l_start+2,bi_next_real,bi_next_imag,res_real,res_imag);
        }*/
        //End iff error check
      }
      //End iff exit
      
      //End iff t
    }
    //End l
  }
  //End recursion step
  int DOES_BACK_CHECK = _TRUE_;
  if(DOES_BACK_CHECK){
    double res_real,res_imag;
    bessel_integral_transform(0,nu_real,nu_imag,t,&res_real,&res_imag);
    double abi_0_abs_sq = abi_real[0]*abi_real[0]+abi_imag[0]*abi_imag[0];
    double lambda_real = (res_real*abi_real[0]+res_imag*abi_imag[0])/abi_0_abs_sq;
    double lambda_imag = (res_imag*abi_real[0]-res_real*abi_imag[0])/abi_0_abs_sq;
    //printf("Found lambda = %.10e+%.10ej \n",lambda_real,lambda_imag);
    
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      temp_real = abi_real[index_l]*lambda_real-abi_imag[index_l]*lambda_imag;
      temp_imag = abi_real[index_l]*lambda_imag+abi_imag[index_l]*lambda_real;
      abi_real[index_l] = temp_real;
      abi_imag[index_l] = temp_imag;
      if(temp_real*temp_real+temp_imag*temp_imag<BESSEL_EPSILON*initial_abs[index_l]){
        //printf("Exiting :: %4d \n",index_l);
        max_t[index_l]=t;
      }
    }
  }
  //clock_t end = clock();
  //printf(" -> Backward Simple took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
  return _SUCCESS_;
}


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
                                            ){
  //clock_t start = clock();
  double bi_real,bi_imag,bi_next_real,bi_next_imag,bi_next_next_real,bi_next_next_imag;
  double temp_real,temp_imag;
  int bessel_recursion_backward_min_l_step = (int)(nu_fraction*(bessel_recursion_backward_min_l_step_high_nu-bessel_recursion_backward_min_l_step_low_nu)+bessel_recursion_backward_min_l_step_low_nu);
  //if(index_coeff>105 && t>0.945){printf("step = %d",bessel_recursion_backward_min_l_step);}
  double min_factor = bessel_recursion_backward_min_l_step*T_SWITCH_FORWARD_RECURSION_SIMPLE;
  int bessel_recursion_l_step = (int)MIN(min_factor/t,bessel_recursion_backward_max_l_step);
  
  int num_recursion_steps = (int)ceil((double)(l_recursion_max+1)/(double)bessel_recursion_l_step);
  int index_recursion_step;
  for(index_recursion_step = 0 ;index_recursion_step < num_recursion_steps;++index_recursion_step){
    int l_start = MAX(l_recursion_max-index_recursion_step*bessel_recursion_l_step,2);
    //printf("Starting at %4d \n",l_start);
    //while(t<max_t[l_start]){
    //  l_start--;
    //}
    //printf("Current rec_step  = %4d , l_max = %4d , l_start = %4d , dl = %4d \n",index_recursion_step,l_max,l_start,(index_recursion_step+1)*bessel_recursion_l_step);
    
    bessel_integral_transform(
                       l_start,
                       nu_real,
                       nu_imag,
                       t,
                       &bi_real,
                       &bi_imag);
    bessel_integral_transform(
                       l_start-1,
                       nu_real,
                       nu_imag,
                       t,
                       &bi_next_real,
                       &bi_next_imag);
    
    abi_real[l_start] = bi_real;
    abi_imag[l_start] = bi_imag;
    abi_real[l_start-1] = bi_next_real;
    abi_imag[l_start-1] = bi_next_imag;
    int index_l;
    for(index_l=l_start;index_l>=(l_recursion_max-(index_recursion_step+1)*bessel_recursion_l_step) && index_l>=2;--index_l){  

      //printf("Allowed, since max_t[%i]=%.10e \n",index_l,max_t[index_l]);
      double l = (double)index_l;

      //printf("Starting with %.10e+%.10ej (l %2d) and %.10e+%.10ej (l-1 %2d) \n",
      //       bi_real,bi_imag,(int)l_start,bi_next_real,bi_next_imag,(int)l_start-1);
      //Explicitly t=1 is important for criterium of 'vanishing' function 
      
      double den_real = (l-2.0+nu_real*0.5);
      double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
      double bi_factor_real = -((1.0+l-nu_real*0.5)*bi_real+nu_imag*bi_imag*0.5);
      double bi_factor_imag = -((1.0+l-nu_real*0.5)*bi_imag-nu_imag*bi_real*0.5);
      double bi_factor_next_real = (1+t*t)/t*(l-0.5)*bi_next_real;
      double bi_factor_next_imag = (1+t*t)/t*(l-0.5)*bi_next_imag;
      //printf("%4d bi facs = %.10e+%.10ej and %.10e+%.10ej \n",index_l,bi_factor_real,bi_factor_imag,bi_factor_next_real,bi_factor_next_imag);
      
      temp_real = (bi_factor_next_real+bi_factor_real);
      temp_imag = (bi_factor_next_imag+bi_factor_imag);
      bi_next_next_real = (temp_real*den_real+temp_imag*nu_imag*0.5)/den_abs_squared;
      bi_next_next_imag = (temp_imag*den_real-temp_real*nu_imag*0.5)/den_abs_squared;
      bi_real = bi_next_real;
      bi_imag = bi_next_imag;
      bi_next_real = bi_next_next_real;
      bi_next_imag = bi_next_next_imag;
      
      abi_real[index_l-2] = bi_next_next_real;
      abi_imag[index_l-2] = bi_next_next_imag;
      /*if(bessel_recursion_error_check>1){
        bessel_integral_transform(l-2.0,nu_real,nu_imag,t,&res_real,&res_imag);
        double err = sqrt((((bi_next_real-res_real)*(bi_next_real-res_real))
              +((bi_next_imag-res_imag)*(bi_next_imag-res_imag)))
              /(res_real*res_real+res_imag*res_imag));
        
        //printf("Err %.10e :: %.10e+%.10ej vs %.10e+%.10ej \n",
        //       err,bi_next_real,bi_next_imag,res_real,res_imag);
        
        class_test(err>bi_allowed_error,
              errmsg,
              "recursion MIDDLE on backward recursion exceeds maximum allowed error at l = %4d, t = %.10e (err = %.10e, seed offset %1d , last seed = %2d) \n Compare wrong %.10e+%.10ei to wanted %.10e+%.10ei",
              (int)l-2,t,err,index_recursion_step,(int)l_start+2,bi_next_real,bi_next_imag,res_real,res_imag);
      }*/
      if(bi_next_real*bi_next_real+bi_next_imag*bi_next_imag<BESSEL_EPSILON*initial_abs[index_l]){ 
        /*max_t[index_l]=t;
        if(index_l==*l_max_cur){
          (*l_max_cur)--;
        }*/
        /*if(bessel_recursion_error_check>0){
          bessel_integral_transform(l-2.0,nu_real,nu_imag,t,&res_real,&res_imag);
          double err = sqrt((((bi_next_real-res_real)*(bi_next_real-res_real))
              +((bi_next_imag-res_imag)*(bi_next_imag-res_imag)))
              /(res_real*res_real+res_imag*res_imag));
          class_test(err>bi_allowed_error,
              errmsg,
              "recursion exit on backward recursion exceeds maximum allowed error at l = %4d, t = %.10e (err = %.10e, seed offset %1d , last seed = %2d) \n Compare wrong %.10e+%.10ei to wanted %.10e+%.10ei",
              (int)l-2,t,err,index_recursion_step,(int)l_start+2,bi_next_real,bi_next_imag,res_real,res_imag);
        }*/
        //End iff error check
      }
      //End iff exit
      
      //End iff t
    }
    //End l
  }
  //End recursion step
  int DOES_BACK_CHECK = _TRUE_;
  if(DOES_BACK_CHECK){
    double res_real,res_imag;
    bessel_integral_transform(0,nu_real,nu_imag,t,&res_real,&res_imag);
    double abi_0_abs_sq = abi_real[0]*abi_real[0]+abi_imag[0]*abi_imag[0];
    double lambda_real = (res_real*abi_real[0]+res_imag*abi_imag[0])/abi_0_abs_sq;
    double lambda_imag = (res_imag*abi_real[0]-res_real*abi_imag[0])/abi_0_abs_sq;
    //printf("Found lambda = %.10e+%.10ej \n",lambda_real,lambda_imag);
    
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      temp_real = abi_real[index_l]*lambda_real-abi_imag[index_l]*lambda_imag;
      temp_imag = abi_real[index_l]*lambda_imag+abi_imag[index_l]*lambda_real;
      abi_real[index_l] = temp_real;
      abi_imag[index_l] = temp_imag;
      if(temp_real*temp_real+temp_imag*temp_imag<BESSEL_EPSILON*initial_abs[index_l]){
        //printf("Exiting :: %4d \n",index_l);
        max_t[index_l]=t;
      }
    }
  }
  //clock_t end = clock();
  //printf(" -> Backward Simple took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
  return _SUCCESS_;
}

int bessel_integral_recursion_forward_simple(int l_max, 
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double nu_fraction,
                                             int bessel_recursion_forward_min_l_step_high_nu,
                                             int bessel_recursion_forward_min_l_step_low_nu,
                                             int bessel_recursion_forward_max_l_step,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             ErrorMsg errmsg
                                             ){
  //double res_real,res_imag;
  //clock_t start = clock();
  double temp_real,temp_imag;
  double bi_real,bi_imag,bi_next_real,bi_next_imag,bi_next_next_real,bi_next_next_imag;
  
  int bessel_recursion_forward_min_l_step = (int)(bessel_recursion_forward_min_l_step_low_nu+(bessel_recursion_forward_min_l_step_high_nu-bessel_recursion_forward_min_l_step_low_nu)*nu_fraction);
  double min_factor = bessel_recursion_forward_min_l_step*sqrt(1-T_SWITCH_FORWARD_RECURSION_SIMPLE);
  int bessel_recursion_l_step = (int)MIN(min_factor/sqrt(1-t),bessel_recursion_forward_max_l_step);

  int num_recursion_steps = (int)ceil((double)(l_max+1)/(double)bessel_recursion_l_step);
  //printf("NUM STEPS = %4d \n",num_recursion_steps);
  int index_recursion_step;
  for(index_recursion_step = 0 ;index_recursion_step < num_recursion_steps;++index_recursion_step){
    int l_start = index_recursion_step*bessel_recursion_l_step;
    
    bessel_integral_transform(
                       l_start,
                       nu_real,
                       nu_imag,
                       t,
                       &bi_real,
                       &bi_imag);
    bessel_integral_transform(
                       l_start+1,
                       nu_real,
                       nu_imag,
                       t,
                       &bi_next_real,
                       &bi_next_imag);
                        
    //printf("Starting with %.10e+%.10ej and %.10e+%.10ej \n",bi_real,bi_imag,bi_next_real,bi_next_imag);
    abi_real[l_start] = bi_real;
    abi_imag[l_start] = bi_imag;
    abi_real[l_start+1] = bi_next_real;
    abi_imag[l_start+1] = bi_next_imag;
    int index_l;
    for(index_l=l_start;index_l<(index_recursion_step+1)*bessel_recursion_l_step;++index_l){  
      if(index_l>l_max){break;}
      /*if(t<max_t[index_l]){
        res_real=0.0;
        res_imag=0.0;
        break;
      }
      else{*/
        abi_real[index_l] = bi_real;
        abi_imag[index_l] = bi_imag;

        if(bi_real*bi_real+bi_imag*bi_imag<BESSEL_EPSILON*initial_abs[index_l]){
          max_t[index_l]=t;
        }
        double l = (double)index_l;
        
        //Explicitly t=1 is important for criterium of 'vanishing' function 
        double den_real = (3.0+l-nu_real*0.5);
        double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
        double bi_factor_real = -((l+nu_real*0.5)*bi_real-nu_imag*bi_imag*0.5);
        double bi_factor_imag = -((l+nu_real*0.5)*bi_imag+nu_imag*bi_real*0.5);
        double bi_factor_next_real = (1+t*t)/t*(l+1.5)*bi_next_real;
        double bi_factor_next_imag = (1+t*t)/t*(l+1.5)*bi_next_imag;
        temp_real = (bi_factor_next_real+bi_factor_real);
        temp_imag = (bi_factor_next_imag+bi_factor_imag);
        bi_next_next_real = (temp_real*den_real-temp_imag*nu_imag*0.5)/den_abs_squared;
        bi_next_next_imag = (temp_imag*den_real+temp_real*nu_imag*0.5)/den_abs_squared;
        
        
        bi_real = bi_next_real;
        bi_imag = bi_next_imag;
        bi_next_real = bi_next_next_real;
        bi_next_imag = bi_next_next_imag;
        
        //End exit check
      //}
      //End iff t
    }
    //End l
  }
  //End recursion step
  //clock_t end = clock();
  //printf(" -> Forward Simple took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
  return _SUCCESS_;
}
/*void hypergeom_short_taylor_series(double a_real,double a_imag,double b_real,double b_imag,double c_real,double c_imag,double z,double* result_real,double* result_imag){
  double C_real = 1.0; double C_imag = 0.0; //< Current term in sum
	double S_real = 1.0; double S_imag = 0.0; //< Current sum accumulation
	double temp_real = C_real;
	double temp_imag = C_imag;
	double ab_real,ab_imag;
	double C_factor;
	double den_abs_squared;
	for (int i=0;i < 5;++i){ //Needs more stopping terms
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		den_abs_squared = (c_real + i)*(c_real + i) + c_imag*c_imag;
		C_factor = z / (i + 1) / den_abs_squared;
		//Result
		C_real = (temp_real*(c_real + i) + temp_imag*c_imag)*C_factor;
		C_imag = (temp_imag*(c_real + i) - temp_real*c_imag)*C_factor;
		S_real += C_real;
		S_imag += C_imag;
	}
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}*/
int bessel_integral_recursion_inverse_self(int l_max,
                                           double nu_real,
                                           double nu_imag,
                                           double t,
                                           double* abi_real,
                                           double* abi_imag,
                                           double* max_t,
                                           double* initial_abs,
                                           ErrorMsg errmsg){
  //clock_t start = clock();
  double res_real,res_imag;
  double pref_first_real,pref_first_imag;
  double pref_scnd_real,pref_scnd_imag;
  double tot_pref_real,tot_pref_imag;
  double temp_real,temp_imag;
  double div_real,div_imag;
  //Recurrent transformed factors of z
  double z = t*t;
  //printf("High t = %.10e ||| ",t);
  double z_fac1 = ((1-z)*(1-z))/((1+z)*(1+z));
  double z_fac2 = 2*(1+z)/(1-z);
  double lnz_fac3 = log(2.0/(1.0+z));
  double l = 0;
  if(nu_imag<NU_IMAG_LIMIT){
    double exp_factor = exp(lnz_fac3*(l+nu_real*0.5));
    double exp_phase = lnz_fac3*nu_imag*0.5;
    double exp_real = exp_factor*cos(exp_phase);
    double exp_imag = exp_factor*sin(exp_phase);
    
    //Calculate total prefactor
    double pref_den_gamma_real = 0.0; double pref_den_gamma_imag = 0.0;
    gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &pref_den_gamma_real, &pref_den_gamma_imag);
    double pref_den_gamma_abs_squared = (pref_den_gamma_real*pref_den_gamma_real+pref_den_gamma_imag*pref_den_gamma_imag);
    tot_pref_real = pow(t,l)*MATH_PI_3_2*(exp_real*pref_den_gamma_real+exp_imag*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
    tot_pref_imag = pow(t,l)*MATH_PI_3_2*(exp_imag*pref_den_gamma_real-exp_real*pref_den_gamma_imag)/pref_den_gamma_abs_squared;
  
    double a_1 = 0.0; double b_1 = 0.0;
    gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
    double a_2 = 0.0; double b_2 = 0.0;
    gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
    double a_3 = 0.0; double b_3 = 0.0;
    gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);
    double abs_3 = a_3*a_3 + b_3*b_3;
    exp_real = (a_1*a_3+b_1*b_3)/abs_3;
    exp_imag = (b_1*a_3-a_1*b_3)/abs_3;
    pref_first_real = a_2*exp_real - b_2*exp_imag;
    pref_first_imag = a_2*exp_imag + b_2*exp_real;
    
    //Calculate prefactor of second hypergeom
    exp_factor = pow(z_fac2, nu_real-2.0);
    double temp = log(z_fac2)*nu_imag;
    double pre_pref_scnd_real = exp_factor*cos(temp);
    double pre_pref_scnd_imag = exp_factor*sin(temp);
    
    double pref_num_gamma_real = 0.0; double pref_num_gamma_imag = 0.0;
    gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &pref_num_gamma_real, &pref_num_gamma_imag);
    
    pref_scnd_real = pre_pref_scnd_real*pref_num_gamma_real - pre_pref_scnd_imag*pref_num_gamma_imag;
    pref_scnd_imag = pre_pref_scnd_real*pref_num_gamma_imag + pre_pref_scnd_imag*pref_num_gamma_real;
  }else{
    //Calculate prefactor including t^l
    double exp_factor = exp(lnz_fac3*(l+nu_real*0.5)+log(t)*l);
    double exp_phase = lnz_fac3*nu_imag*0.5;
    double exp_real = exp_factor*cos(exp_phase);
    double exp_imag = exp_factor*sin(exp_phase);
    
    tot_pref_real = MATH_PI_3_2*exp_real;
    tot_pref_imag = MATH_PI_3_2*exp_imag;
    
    //Calculate gamma factors (large nu -> lnGamma)
    double a_0 = 0.0; double b_0 = 0.0;
    ln_gamma_complex_2(1.5-nu_real*0.5, -nu_imag*0.5, &a_0, &b_0);
    double a_1 = 0.0; double b_1 = 0.0;
    ln_gamma_complex_2(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
    double a_2 = 0.0; double b_2 = 0.0;
    ln_gamma_complex_2(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
    double a_3 = 0.0; double b_3 = 0.0;
    ln_gamma_complex_2(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);

    exp_factor = exp(a_1+a_2 - a_3 -a_0);
    exp_real = exp_factor*cos(b_1+b_2 - b_3-b_0);
    exp_imag = exp_factor*sin(b_1+b_2 - b_3-b_0);
    
    //Calculate first prefactor
    pref_first_real = exp_real;
    pref_first_imag = exp_imag;
    
    //Calculate prefactor of second hypergeom
    exp_factor = pow(z_fac2, nu_real-2.0);
    double temp = log(z_fac2)*nu_imag;
    double pre_pref_scnd_real = exp_factor*cos(temp);
    double pre_pref_scnd_imag = exp_factor*sin(temp);
    
    double c_0 = 0.0; double d_0 = 0.0;
    ln_gamma_complex_2(nu_real*0.5 - 1, nu_imag*0.5, &c_0, &d_0);
    exp_factor = exp(c_0-a_0);
    exp_real = exp_factor*cos(d_0-b_0);
    exp_imag = exp_factor*sin(d_0-b_0);
    
    //Calculate second prefactor
    double pref_num_gamma_real = exp_real;
    double pref_num_gamma_imag = exp_imag;
    pref_scnd_real = pre_pref_scnd_real*pref_num_gamma_real - pre_pref_scnd_imag*pref_num_gamma_imag;
    pref_scnd_imag = pre_pref_scnd_real*pref_num_gamma_imag + pre_pref_scnd_imag*pref_num_gamma_real;
  }
  int index_l;
  for(index_l=0;index_l<=l_max;++index_l){
    double l = (double)index_l;
    //Calculate first hypergeom  
    double first_hypergeom_real = 0.0; double first_hypergeom_imag = 0.0;
    hypergeom(l*0.5 + nu_real*0.25, nu_imag*0.25, nu_real*0.25 + l*0.5 + 0.5, nu_imag*0.25, nu_real*0.5, nu_imag*0.5, z_fac1, &first_hypergeom_real, &first_hypergeom_imag);
    
    //Calculate second hypergeom
    double scnd_hypergeom_real = 0.0; double scnd_hypergeom_imag = 0.0;
    hypergeom(l*0.5 +1.5 - nu_real*0.25, -nu_imag*0.25, 1.0 + l*0.5 - nu_real*0.25, -nu_imag*0.25, 2 - nu_real*0.5, -nu_imag*0.5, z_fac1, &scnd_hypergeom_real, &scnd_hypergeom_imag);
    
    //Calculate sum of hypergeoms with their prefactors
    double total_bracket_real = (first_hypergeom_real*pref_first_real - first_hypergeom_imag*pref_first_imag) + (scnd_hypergeom_real*pref_scnd_real - scnd_hypergeom_imag*pref_scnd_imag);
    double total_bracket_imag = (first_hypergeom_real*pref_first_imag + first_hypergeom_imag*pref_first_real) + (scnd_hypergeom_real*pref_scnd_imag + scnd_hypergeom_imag*pref_scnd_real);
    
    //Calculate total
    res_real = total_bracket_real*tot_pref_real - total_bracket_imag*tot_pref_imag;
    res_imag = total_bracket_real*tot_pref_imag + total_bracket_imag*tot_pref_real;
    
    
    tot_pref_real*=(2./(1.+z))*t;
    tot_pref_imag*=(2./(1.+z))*t;
    
    div_real = l+1-nu_real/2. +1;
    div_imag = -nu_imag/2.;
    double div_abs_sq = div_real*div_real+div_imag*div_imag;
    double factor_real = ((l-1+nu_real*0.5 +1)*div_real+div_imag*nu_imag*0.5)/div_abs_sq;
    double factor_imag = (nu_imag*0.5*div_real-(1+   l-1+nu_real*0.5)*div_imag)/div_abs_sq;
     
    temp_real = factor_real*pref_first_real-factor_imag*pref_first_imag;
    temp_imag = factor_real*pref_first_imag+factor_imag*pref_first_real;
    
    pref_first_real = temp_real;
    pref_first_imag = temp_imag;
    //printf("HELLO ! \n");
    //bessel_integral_transform(l,nu_real,nu_imag,t,&res_real,&res_imag);
    if(res_real*res_real+res_imag*res_imag<BESSEL_EPSILON*initial_abs[index_l]){
      max_t[index_l]=t;
    }
    abi_real[index_l] = res_real;
    abi_imag[index_l] = res_imag;
    //End exit check
  }
  //End l
  //clock_t end = clock();
  //printf(" -> Inverse Self took %.10e seconds \n",((double)(end-start))/CLOCKS_PER_SEC);
  return _SUCCESS_;
}
