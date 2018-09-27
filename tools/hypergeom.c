/** @file functions.c Documented special functions module.
 *
 * Nils Schoenberg, 12.10.2017
 *
 * This module computes complex gamma functions, the gaussian hypergeometric function (Hypergeometric2F1) and various derived functions (bessel integrals)
 * Additionally, it uses various types of recursions.
 */
#include "common.h"
#include "hypergeom.h"
#include <math.h>



/* Mathematical constants used in this file */

//Variations of PI
#define MATH_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930382 /** < PI , 200 digits*/
#define MATH_PI_3_2 5.5683279968317078452848179821188357020136243902832439107536758188297455336477957022121776873847084940970621035598961308638949212663157851705967389211068321811703451381324726069689321738560369691961862 /** < PI^(3/2) , 200 digits*/
#define MATH_2_PI 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639 /** < 2 PI , 200 digits*/
#define MATH_PI_HALF 1.570796326794896619231321691639751442098584699687552910487472296153908203143104499314017412671058533991074043256641153323546922304775291115862679704064240558725142051350969260552779822311474477465191/** < PI/2 , 200 digits */
#define SQRT_2_PI 2.5066282746310005024157652848110452530069867406099383166299235763422936546078419749465958383780572661160099726652038796448663236181267361809578556559661479319134354804581237311373278043130199388026472 /** < sqrt(2*PI) , 200 digits */

//Logarithms of PI
#define ln2PI 1.8378770664093454835606594728112352797227949472755668256343030809655313918545207953894865972719083952440112932492686748927337257636815871443117518304453627872071214850947173380927918119827616112603265 /** < ln(2*PI) , 200 digits */
#define lnPI_HALF_05 0.22579135264472743236309761494744107178589733927752815869647153098937207395756568208887997163953551008000416560406365171268134264608266301512320516358760543167283317171898723385608886736748913215903988 /** < ln(PI/2)/2 , 200 digits */

//Logarithms
#define LOG_2 0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102057068573368552023575813055703267075163507596193072757082837143519030703862389167347112335 /** < ln(2) , 200 digits*/



/* Constants definition for the Lanczos Approximation */
const int NUM_GAMMA_COEFF = 8;
const double GAMMA_COEFF[] = { 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };



/* Precision parameters */
double BESSEL_INTEGRAL_T_SWITCH = 0.95;/** < switching value of bessel integrals, low t safe to high t, Default : 0.95*/
double BESSEL_INTEGRAL_T_SWITCH_SAFE = 0.9;/** < switching value of bessel integrals, low t to low t safe, Default : 0.9*/
double TAYLOR_HYPERGEOM_ACCURACY = 1e-10; /** < numerical accuracy */
const int HYPERGEOM_MAX_ITER = 100000; /** < Maximum number of iterations in taylor sums */
double ABS_LIMIT = 1e100; /** < Maximum acceptable absolute magnitude for a number, whose square should (easily) fit double precision */
double NU_IMAG_LIMIT = 150; /** < Maximum value of nu_imag such that cosh(nu_imag) or gamma(nu_imag) doesn't overflow */
double BESSEL_INTEGRAL_MAX_L = 80.0;/** < Maximum acceptable l for Gamma(l+...) to be still finite */
double MAX_SINE_IMAG = 100.0; /** < Maximum value that is acceptable for cosh(z) doesn't overflow */



/* Function declarations */
void hypergeom(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag);



/* Function implementations */
/**
 * This small helper function computes the logarithm of the sine for large imaginary part of the argument (y)
 * Basically we ignore the small exponent of exp(i*(x+iy))=exp(-y+ix) and exp(-i*(x+iy))=exp(y-ix)
 * Depending on the sign of z_imag we can ignore one of these
 * 
 * ! Assumes abs(z_imag) big
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
 * This small helper function computes the logarithm of the cosine for large imaginary part of the argument (y)
 * Basically we ignore the small exponent of exp(i*(x+iy))=exp(-y+ix) and exp(-i*(x+iy))=exp(y-ix)
 * Depending on the sign of z_imag we can ignore one of these
 * 
 * ! Assumes abs(z_imag) big
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
 * This function computes the real Gamma function, which can be found in the system library
 * */
double gamma_real(double z){
	return tgamma(z);
}
/**
 * This function computes the logarithm of the real Gamma function, which can be found in the system library
 * */
double ln_gamma_real(double z){
	return lgamma(z);
}
/**
 * This function computes the complex Gamma function 
 *  using the Lanczos Approximation to order of NUM_GAMMA_COEFF,
 *  using the GAMMA_COEFF array as coefficients
 * The analytical continuation to Re(z)<0.5 is done using the reflection formula of the gamma function
 * Turns out to be a bit faster than the sinh(1/z) approximation
 * 
 * The formula is 
 * Gamma(z+1) = sqrt(2pi) t^(z+1/2) e^(-t) sum c_n ((z(z-1)...(z-n+1))/((z+1)(z+2)...(z+n))
 *   where t = z+N-1/2 in our case
 * See also https://en.wikipedia.org/wiki/Lanczos_approximation
 * 
 * 
 * ! Assumes MATH_PI*z_imag to be an overflow-safe exponent
 * */
void gamma_complex(double z_real,double z_imag,double* res_real,double* res_imag){
	if (z_real > 0.5){
		z_real -= 1;
		//Initialize Lanczos sum
		double sum_real = 0.99999999999980993;
		double sum_imag = 0.0;
		double temp_abs_squared; double temp_real;
    int i;
    //Sum over the coefficients
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * temp_real / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t=z+N-1/2 and log(t)
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
    //Sum over the coefficients
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * (z_real + i + 1) / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t=z+N-1/2 and log(t)
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
 * This function computes the logarithm of the complex Gamma function 
 *  using the Lanczos Approximation to order of NUM_GAMMA_COEFF,
 *  using the GAMMA_COEFF array as coefficients
 * The analytical continuation to Re(z)<0.5 is done using the reflection formula of the gamma function
 * Turns out to be a bit faster than the sinh(1/z) approximation
 * 
 * The formula is 
 * Gamma(z+1) = sqrt(2pi) t^(z+1/2) e^(-t) sum c_n ((z(z-1)...(z-n+1))/((z+1)(z+2)...(z+n))
 *   where t = z+N-1/2 in our case
 * See also https://en.wikipedia.org/wiki/Lanczos_approximation
 * */
void ln_gamma_complex(double z_real,double z_imag,double* res_real,double* res_imag){
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
		//Calculate t=z+N-1/2, log(t) and log(sum)
		double t_real = z_real + NUM_GAMMA_COEFF - 0.5;
		double log_sum_real = 0.5*log(sum_real*sum_real+sum_imag*sum_imag);
		double log_sum_imag = atan2(sum_imag,sum_real);
		double log_t_real = 0.5*log(t_real*t_real+z_imag*z_imag);
		double log_t_imag = atan2(z_imag,t_real);
		//Return result
		*res_real = 0.5*ln2PI+log_sum_real-t_real+log_t_real*(0.5+z_real)-log_t_imag*z_imag;
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
    //Sum over coefficents
		for (i = 0; i < NUM_GAMMA_COEFF; ++i){
			temp_real = z_real+i+1;
			temp_abs_squared = temp_real * temp_real + z_imag*z_imag;
			sum_real += GAMMA_COEFF[i] * temp_real / temp_abs_squared;
			sum_imag -= GAMMA_COEFF[i] * (z_imag) / temp_abs_squared;
		}
		//Calculate t=z+N-1/2, log(t) and log(sum)
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
 * This function calculates the Gaussian hpyergeometric function 
 *  for a,b complex, c,z real as a simple Taylor series
 * 
 * ! Assumes that no overflow occurs during the hypergeometric summation
 * */
void hypergeom_series_sum_ab(double a_real, double a_imag, double b_real, double b_imag, double c, double z, double* result_real, double* result_imag){
  //Initial variable setup
	double C_real = 1.0; double C_imag = 0.0; //< Current term in sum
	double S_real = 1.0; double S_imag = 0.0; //< Current sum accumulation
	int i = 0;
	double temp_real = C_real; double temp_imag = C_imag; //< temporary spaceholder
	double ab_real,ab_imag;
	double C_factor;
	//Stop when convergence or maximum iterations reached
	while ((C_real*C_real + C_imag*C_imag) > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY*(S_real*S_real + S_imag*S_imag)  && i <= HYPERGEOM_MAX_ITER){
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		C_factor = z / (i + 1) / (c + i);
		//Add term to sum
		C_real = temp_real*C_factor;
		C_imag = temp_imag*C_factor;
		S_real += C_real;
		S_imag += C_imag;
		i += 1;
	}
  //Return result
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the Gaussian hpyergeometric function 
 *  for a,b,c complex, z real as a simple Taylor series
 * 
 * ! Assumes that no overflow occurs during the hypergeometric summation
 * */
void hypergeom_series_sum_abc(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag){
  //Initial variable setup
	double C_real = 1.0; double C_imag = 0.0; //< Current term in sum
	double S_real = 1.0; double S_imag = 0.0; //< Current sum accumulation
	int i = 0;
	double temp_real = C_real; double temp_imag = C_imag; //< temporary spaceholder
	double ab_real,ab_imag;
	double C_factor;
	double den_abs_squared;
	//Stop when convergence or maximum iterations reached
	while ((C_real*C_real + C_imag*C_imag) > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY *(S_real*S_real + S_imag*S_imag) && i <= HYPERGEOM_MAX_ITER){
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		den_abs_squared = (c_real + i)*(c_real + i) + c_imag*c_imag;
		C_factor = z / (i + 1) / den_abs_squared;
		//Add term to sum
		C_real = (temp_real*(c_real + i) + temp_imag*c_imag)*C_factor;
		C_imag = (temp_imag*(c_real + i) - temp_real*c_imag)*C_factor;
		S_real += C_real;
		S_imag += C_imag;
		i += 1;
	}
  //Return result
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the Gaussian hpyergeometric function 
 *  for a,b complex, c,z real as a simple Taylor series
 * 
 * Overflow-safe version, stores number of overflows of size 
 *  ABS_LIMIT in the overflows variable.
 * */
void hypergeom_series_sum_ab_safe(double a_real, double a_imag, double b_real, double b_imag, double c, double z, double* result_real, double* result_imag, int* overflows){
	double C_real = 1.0; double C_imag = 0.0; //< Current term in sum
	double S_real = 1.0; double S_imag = 0.0; //< Current sum accumulation
	int i = 0;
  *overflows = 0;
	double temp_real = C_real; double temp_imag = C_imag; //< temporary spaceholder
	double ab_real,ab_imag;
	double C_factor;
  double C_abs = 1.0,S_abs=1.0;
	//Stop when convergence or maximum iterations reached
	while (C_abs > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY*S_abs && i <= HYPERGEOM_MAX_ITER){
		//Calculate numerator
    ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
    ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
    temp_real = C_real*ab_real - C_imag*ab_imag;
    temp_imag = C_imag*ab_real + C_real*ab_imag;
    //Calculate denominator
    C_factor = z / (i + 1) / (c + i);
    C_real = temp_real*C_factor;
    C_imag = temp_imag*C_factor;
    //Calculate new absolute values and check for divergence
    C_abs = C_real*C_real+C_imag*C_imag;
    S_abs = S_real*S_real+S_imag*S_imag;
    if(C_abs>ABS_LIMIT*ABS_LIMIT || S_abs > ABS_LIMIT*ABS_LIMIT){
      //If summation would overflow, instead divide system through large number,
      // and continue with overflows increased by one instead
      C_real/=ABS_LIMIT;
      C_imag/=ABS_LIMIT;
      S_real/=ABS_LIMIT;
      S_imag/=ABS_LIMIT;
      (*overflows)+=1;
    }
    //Add new term to sum
    S_real += C_real;
    S_imag += C_imag;
    i+=1;
	}
  //Return result
  *result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the Gaussian hpyergeometric function 
 *  for a,b,c complex, z real as a simple Taylor series
 * 
 * Overflow-safe version, stores number of overflows of size 
 *  ABS_LIMIT in the overflows variable.
 * */
void hypergeom_series_sum_abc_safe(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag,int* overflows){
	double C_real = 1.0; double C_imag = 0.0; // Current term in sum
	double S_real = 1.0; double S_imag = 0.0; // Current sum accumulation*
	int i = 0;
  *overflows = 0;
	double temp_real = C_real; double temp_imag = C_imag; //< temporary spaceholder
	double ab_real,ab_imag;
	double C_factor;
  double C_abs = 1.0,S_abs=1.0;
	double den_abs_squared;
	//Stop when convergence or maximum iterations reached
	while (C_abs > TAYLOR_HYPERGEOM_ACCURACY*TAYLOR_HYPERGEOM_ACCURACY *S_abs && i <= HYPERGEOM_MAX_ITER){
		//Calculate numerator
		ab_real = (a_real+i)*(b_real+i)-a_imag*b_imag;
		ab_imag = (a_real+i)*b_imag+a_imag*(b_real+i);
		temp_real = C_real*ab_real - C_imag*ab_imag;
		temp_imag = C_imag*ab_real + C_real*ab_imag;
		//Calculate denominator
		den_abs_squared = (c_real + i)*(c_real + i) + c_imag*c_imag;
		C_factor = z / (i + 1) / den_abs_squared;
		C_real = (temp_real*(c_real + i) + temp_imag*c_imag)*C_factor;
		C_imag = (temp_imag*(c_real + i) - temp_real*c_imag)*C_factor;
    //Calculate new absolute values and check for divergence
    C_abs = C_real*C_real+C_imag*C_imag;
    S_abs = S_real*S_real+S_imag*S_imag;
    if(C_abs>ABS_LIMIT*ABS_LIMIT || S_abs > ABS_LIMIT*ABS_LIMIT){
      //If summation would overflow, instead divide system through large number,
      // and continue with overflows increased by one instead
      C_real/=ABS_LIMIT;
      C_imag/=ABS_LIMIT;
      S_real/=ABS_LIMIT;
      S_imag/=ABS_LIMIT;
      (*overflows)+=1;
    }
    //Add new term to sum
		S_real += C_real;
		S_imag += C_imag;
		i += 1;
	}
  //Return result
	*result_real = S_real;
	*result_imag = S_imag;
	return;
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * It selects from previously defined functions, and handles special cases
 *
 * This function is not guaranteed to give a reasonable result, especially for z->1.
 * The handling of the convergence is instead done by the bessel_integral functions!
 * 
 * ! Assumes z < 1
 * ! Overflow can occur
 * */
void hypergeom(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag){
	//Should the border value of z at 0 be used, we know the value of the hypergeometric function immediately
	if(z==0.0){
		*result_real = 1; *result_imag=0; return;
	}
	//In the case of c_imag ==0.0, we should use the simplified version to save calculation time
	if(c_imag==0.0){
		hypergeom_series_sum_ab(a_real,a_imag,b_real,b_imag,c_real,z,result_real,result_imag);
		return;
	}
	//In the general case, we use the general version of the taylor approximation
	else{
		hypergeom_series_sum_abc(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag);
		return;
	}
}
/**
 * This function calculates the hpyergeometric gaussian function for a,b,c complex,z real
 * It selects from previously defined functions, and handles special cases, and is overflow-safe
 *
 * This function is not guaranteed to give a reasonable result, especially for z->1.
 * The handling of the convergence is instead done by the bessel_integral functions!
 * 
 * ! Assumes z < 1
 * */
void hypergeom_safe(double a_real, double a_imag, double b_real, double b_imag, double c_real, double c_imag, double z, double* result_real, double* result_imag,int* overflows){
  //Should the border value of z at 0 be used, we know the value of the hypergeometric function immediately
	if(z==0.0){
		*result_real = 1; *result_imag=0; *overflows=0; return;
	}
	//In the case of c_imag ==0.0, we should use the simplified version to save calculation time
	if(c_imag==0.0){
		hypergeom_series_sum_ab_safe(a_real,a_imag,b_real,b_imag,c_real,z,result_real,result_imag,overflows);
		return;
	}
	//In the general case, we use the general version of the taylor approximation
	else{
		hypergeom_series_sum_abc_safe(a_real,a_imag,b_real,b_imag,c_real,c_imag,z,result_real,result_imag,overflows);
		return;
	}
}
/**
 * This function calculates the analytical limit 
 *  for the bessel_integral for t->1 (z=t²->1 as well)
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
		ln_gamma_complex(l + nu_real*0.5, nu_imag*0.5, &first_ln_gamma_real, &first_ln_gamma_imag);
		double scnd_ln_gamma_real = 0.0; double scnd_ln_gamma_imag = 0.0;
		ln_gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &scnd_ln_gamma_real, &scnd_ln_gamma_imag);
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
		ln_gamma_complex(1-nu_real*0.5, -nu_imag*0.5, &first_gamma_real, &first_gamma_imag);
		double scnd_gamma_real = 0.0; double scnd_gamma_imag = 0.0;
		ln_gamma_complex(1.5- nu_real*0.5, -nu_imag*0.5, &scnd_gamma_real, &scnd_gamma_imag);

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
 * This function calculates the integral over the bessel functions at low t slightly differently, using a single transformation
 * This transformation reduces the parameters l and nu (divides them by 2), but increases z (roughly factor of 4)
 *
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1 is not true
 * Here, the CALLING function will have care about always keeping t<<1 .
 *
 * This function is mostly applicable in the low t regime
 * 
 * ! Assumes t << 1
 * ! Can overflow
 * */
void bessel_integral_lowt_transform(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Common parameters
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
		ln_gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
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
		ln_gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
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
/**
 * This function calculates the integral over the bessel functions at low t slightly differently, using a single transformation
 * This transformation reduces the parameters l and nu (divides them by 2), but increases z (roughly factor of 4)
 * 
 * It is overflow-safe 
 * 
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1 is not true
 * Here, the CALLING function will have care about always keeping t<<1 .
 *
 * This function is mostly applicable in the low t regime
 * 
 * ! Assumes t << 1
 * */
void bessel_integral_lowt_transform_safe_flag(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag,short* overflow_flag){
	//Common parameters
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
		ln_gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
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
		ln_gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
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
 * This function calculates the integral over the bessel functions at low t slightly differently, using a single transformation
 * This transformation reduces the parameters l and nu (divides them by 2), but increases z (roughly factor of 4)
 * 
 * It is overflow-safe 
 * 
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1 is not true
 * Here, the CALLING function will have care about always keeping t<<1 .
 *
 * This function is mostly applicable in the low t regime
 * 
 * ! Assumes t << 1
 * */
void bessel_integral_lowt_transform_safe(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	short overflow_flag = _FALSE_;
  //Does safe evaluation, but ignores the flag
  bessel_integral_lowt_transform_safe_flag(l,nu_real,nu_imag,t,res_real,res_imag,&overflow_flag);
  return;
}
/**
 * This function calculates the integral over the bessel functions for high t, but in a slightly different way
 *
 * VERY IMPORTANT: The function will give INCORRECT results due to numerical errors for values where t<<1
 * Special care has to be taken because of that, which is explained further in the paper.
 * Here, the CALLING function will have care about always keeping (1-t*t)<<1 .
 *
 * This function is mostly applicable in the high t regime
 * 
 * ! Assumes (1-t*t)<<1
 * */
void bessel_integral_hight_transform(double l, double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
	//Common factors
	double z = t*t;
	double z_fac1 = ((1-z)*(1-z))/((1+z)*(1+z));
	double z_fac2 = 2*(1+z)/(1-z);
	double lnz_fac3 = log(2.0/(1.0+z));
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

			//Calculate gamma factors for large l
			double a_1 = 0.0; double b_1 = 0.0;
			ln_gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
			double a_2 = 0.0; double b_2 = 0.0;
			gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
			double a_3 = 0.0; double b_3 = 0.0;
			ln_gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);

			exp_factor = exp(a_1 - a_3);
			exp_real = exp_factor*cos(b_1 - b_3);
			exp_imag = exp_factor*sin(b_1 - b_3);

			//Calculate first prefactor
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
		ln_gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &a_0, &b_0);
		double a_1 = 0.0; double b_1 = 0.0;
		ln_gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
		double a_2 = 0.0; double b_2 = 0.0;
		ln_gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
		double a_3 = 0.0; double b_3 = 0.0;
		ln_gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);

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
		ln_gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &c_0, &d_0);
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
/**
 * The I_l(nu,t) for l = 0 is known analytically. The below implements the corresponding formula.
 *  There is a small expansion required for t<<1 due to numerical divergences of the full formula
 * */
#define N_BESSEL_L0_COEFFS 16
const double bessel_integral_l0_coeffs[N_BESSEL_L0_COEFFS]={4.,2./3.,1./30., 1.0/1260.,1./90720.,1./9979200.,1./1556755200.,1./326918592000.,1./88921857024000.,1./30411275102208000.,1./12772735542927360000.,1./6463004184721244160000.,1./3877802510832746496000000.,1./2722217362604588040192000000.,1./2210440498434925488635904000000.,1./2055709663544480704431390720000000.};
void bessel_integral_l0(double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
  //Treat the analytic limit at t=1
  if( t == 1 ){
    bessel_analytic_limit_atz1(0,nu_real,nu_imag,res_real,res_imag);
    return;
  }
  /**
   * Calculate the result using
   *  coeff[i]*(2*pi*cos(nu*pi/2)/2*gamma(nu-2)*(nu-2)/2*gamma(nu-2+2*i)/gamma(nu-2) * (2*i+nu)
   * 
   * The coefficients decline even after being multiplied with large gamma functions for t->1
   * As such, the series would theoretically converge even towards t->1,
   *  but there only very VERY slowly. We expect in that regime the (simple) formula below to hold better anyway
   *
   * Thus we here explicitly use this formula only for t<T_MIN_TAYLOR
   * 
   * Description:
   * 
   * Use taylor series expansion for small t to avoid numerical cancelations
   *  t^20 is already at least of relative accuracy 1e-20, 
   *  which is precise enough for double precision variables
   * However, due to numerical roundoffs, it does not quite reach as low
   *  It is still good enough in that case, though
   * */
  else if(t<T_MIN_TAYLOR){
    double pref_real,pref_imag;
    //Overflow-safe version of 2*pi*e^(i*pi*nu/2)*gamma(nu-2)
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
      ln_gamma_complex(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    //Getting the total prefactor
    double tot_pref_real = (nu_real*0.5-1.0)*pref_real-nu_imag*pref_imag*0.5;
    double tot_pref_imag = (nu_real*0.5-1.0)*pref_imag+nu_imag*pref_real*0.5;

    //Setup coefficient summation
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
    //Summing the coefficient formula
    for(i=0;i<N_BESSEL_L0_COEFFS;++i){
      nu_pref_real = prod_nu_fac_real;
      nu_pref_imag = prod_nu_fac_imag;
      cur_real = -bessel_integral_l0_coeffs[i]*cur_t*(nu_pref_real*tot_pref_real-nu_pref_imag*tot_pref_imag);
      cur_imag = -bessel_integral_l0_coeffs[i]*cur_t*(nu_pref_real*tot_pref_imag+nu_pref_imag*tot_pref_real);

      cur_t *= t*t;
      //Formula for next coefficient as a function of previous coefficient
      cur_nu_fac_real = ((nu_real+2.0*i-1.0)*(nu_real+2.0*i)-nu_imag*nu_imag);
      cur_nu_fac_imag = nu_imag*(2.0*nu_real+4.0*i-1.0);
      temp_real = cur_nu_fac_real*prod_nu_fac_real-cur_nu_fac_imag*prod_nu_fac_imag;
      temp_imag = cur_nu_fac_real*prod_nu_fac_imag+cur_nu_fac_imag*prod_nu_fac_real;
      //Add new terms to sum
      prod_nu_fac_real = temp_real;
      prod_nu_fac_imag = temp_imag;
      sum_real += cur_real;
      sum_imag += cur_imag;
    }
    //Get result
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
  /**
   * Calculate the result using the full, but simple formula
   * The simple formula is only applicable, when NOT t<<1
   * */
  else{
    //Common variables
    double log1pt = log(1+t);
    double log1mt = log(1-t);
    double pref_real,pref_imag;
    //Overflow-safe version of 2*pi*e^(i*pi*nu/2)*gamma(nu-2)
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
      ln_gamma_complex(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    //Calculate (1-t)^(2-nu) and (1+t)^(2-nu)
    double t_fac_first = exp(log1pt*(2.0-nu_real));
    double t_fac_second = -exp(log1mt*(2.0-nu_real));
    double t_fac_real = t_fac_first*cos(log1pt*nu_imag)+t_fac_second*cos(log1mt*nu_imag);
    double t_fac_imag = -t_fac_first*sin(log1pt*nu_imag)-t_fac_second*sin(log1mt*nu_imag);
    //Get result
    *res_real = (t_fac_real*pref_real-t_fac_imag*pref_imag)/t;
    *res_imag = (t_fac_real*pref_imag+t_fac_imag*pref_real)/t;
    return;
  }
}
/**
 * The I_l(nu,t) for l = 1 is known analytically. The below implements the corresponding formula.
 *  There is a small expansion required for t<<1 due to numerical divergences of the full formula
 * */
#define N_BESSEL_L1_COEFFS 11
const double bessel_integral_l1_coeffs[N_BESSEL_L1_COEFFS]={4./3.,2./15.,1./210., 1.0/11340.,1./997920.,1./129729600.,1./23351328000.,1./5557616064000.,1./1689515283456000.,1./638636777146368000.,1./293772917487329280000.};
//Interestingly, the above coefficients are related to the l0 coefficients by 1/(2n+1) (with n=index_in_array_starting_from_0+1)
void bessel_integral_l1(double nu_real, double nu_imag, double t, double* res_real, double* res_imag){
  //Treat analytically the t=1 case
  if( t == 1 ){
    bessel_analytic_limit_atz1(1,nu_real,nu_imag,res_real,res_imag);
    return;
  }
  /**
   * Calculate the result using
   *  coeff[i]*(2*pi*cos(nu*pi/2)/2*gamma(nu-2)*(nu-2)/2*gamma(nu-2+2*i)/gamma(nu-2) * (2*i+nu)
   * 
   * The coefficients decline even after being multiplied with large gamma functions for t->1
   * As such, the series would theoretically converge even towards t->1,
   *  but there only very VERY slowly. We expect in that regime the (simple) formula below to hold better anyway
   *
   * Thus we here explicitly use this formula only for t<T_MIN_TAYLOR
   * 
   * Description:
   * 
   * Use taylor series expansion for small t to avoid numerical cancelations
   *  t^21 is already at least of relative accuracy 1e-21, 
   *  which is precise enough for double precision variables
   * However, due to numerical roundoffs, it does not quite reach as low
   *  It is still good enough in that case, though
   * */
  else if(t<T_MIN_TAYLOR){
    double pref_real,pref_imag;
    //Overflow-safe version of 2*pi*e^(i*pi*nu/2)*gamma(nu-2)
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
      ln_gamma_complex(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    //Setup coefficient summation
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
    //Summing the coefficient formula
    for(i=0;i<N_BESSEL_L1_COEFFS;++i){
      nu_pref_real = (nu_real+2.0*i)*prod_nu_fac_real-nu_imag*prod_nu_fac_imag;
      nu_pref_imag = (nu_real+2.0*i)*prod_nu_fac_imag+nu_imag*prod_nu_fac_real;
      cur_real = -bessel_integral_l1_coeffs[i]*cur_t*(nu_pref_real*tot_pref_real-nu_pref_imag*tot_pref_imag);
      cur_imag = -bessel_integral_l1_coeffs[i]*cur_t*(nu_pref_real*tot_pref_imag+nu_pref_imag*tot_pref_real);

      cur_t *= t*t;
      //Formula for next coefficient as a function of previous coefficient
      cur_nu_fac_real = ((nu_real+2.0*i-1.0)*(nu_real+2.0*i)-nu_imag*nu_imag);
      cur_nu_fac_imag = nu_imag*(2.0*nu_real+4.0*i-1.0);
      temp_real = cur_nu_fac_real*prod_nu_fac_real-cur_nu_fac_imag*prod_nu_fac_imag;
      temp_imag = cur_nu_fac_real*prod_nu_fac_imag+cur_nu_fac_imag*prod_nu_fac_real;
      //Add new terms to the sum
      prod_nu_fac_real = temp_real;
      prod_nu_fac_imag = temp_imag;
      sum_real += cur_real;
      sum_imag += cur_imag;
    }
    //Get result
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
  /**
   * Calculate the result using the full, but simple formula
   * The simple formula is only applicable, when NOT t<<1
   * */
  else{
    //Common variables
    double log1pt = log(1+t);
    double log1mt = log(1-t);
    double pref_real,pref_imag;
    //Overflow-safe version of 2*pi*e^(i*pi*nu/2)*gamma(nu-2)
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
      ln_gamma_complex(nu_real-2.0,nu_imag,&ln_gam_real,&ln_gam_imag);
      exp_factor = exp(ln_gam_real+ln_cos_real);
      pref_real = MATH_2_PI*exp_factor*cos(ln_gam_imag+ln_cos_imag);
      pref_imag = MATH_2_PI*exp_factor*sin(ln_gam_imag+ln_cos_imag);
    }
    //Now for the more complicated version of the analytical result, including a 1/(4-nu) prefactor
    double nu_fac_abs_squared = (4.0-nu_real)*(4.0-nu_real)+nu_imag*nu_imag;
    double tot_pref_real = ((4.0-nu_real)*pref_real-nu_imag*pref_imag)/nu_fac_abs_squared;
    double tot_pref_imag = ((4.0-nu_real)*pref_imag+nu_imag*pref_real)/nu_fac_abs_squared;

    //And (1-t)^(2-nu),(1+t)^(2-nu) and (1-t)^2 + nu*t and (1+t)^2 - nu*t
    double t_fac_first = exp(log1pt*(2.0-nu_real));
    double t_fac_second = -exp(log1mt*(2.0-nu_real));
    double t_fac_bracket_first_real = (1-t)*(1-t)+nu_real*t;
    double t_fac_bracket_second_real = (1+t)*(1+t)-nu_real*t;
    double first_real = t_fac_first*(cos(log1pt*nu_imag)*t_fac_bracket_first_real+sin(log1pt*nu_imag)*nu_imag*t);
    double first_imag = t_fac_first*(-sin(log1pt*nu_imag)*t_fac_bracket_first_real+cos(log1pt*nu_imag)*nu_imag*t);
    double second_real = t_fac_second*(cos(log1mt*nu_imag)*t_fac_bracket_second_real-sin(log1mt*nu_imag)*nu_imag*t);
    double second_imag = t_fac_second*(-sin(log1mt*nu_imag)*t_fac_bracket_second_real-cos(log1mt*nu_imag)*nu_imag*t);
    
    //Get result
    *res_real = ((first_real+second_real)*tot_pref_real-(first_imag+second_imag)*tot_pref_imag)/t/t;
    *res_imag = ((first_real+second_real)*tot_pref_imag+(first_imag+second_imag)*tot_pref_real)/t/t;
    return;
  }
}

/**
 * The I_l(nu,t) can still be Taylor-expanded for small t<<1,
 *  corresponding to the hypergeometric taylor sum for t^2 << 1 being 
 *  stopped at some coefficient N_BESSEL_LTAYLOR_COEFFS
 * 
 * Tradeoff: Speed vs accuracy in the choice of this maximum coefficient
 * */
#define N_BESSEL_LTAYLOR_COEFFS 40
void bessel_taylor_for_small_t(double l,double nu_real,double nu_imag,double t, double* res_real,double* res_imag){
  //Bessel integral prefactor
  double pref_t = pow(t,l);
  double pref_real,pref_imag;
  //Evalutae here the prefactor in an overflow-safe way
  if(nu_imag>NU_IMAG_LIMIT){
		double a_0 = 0.0; double  b_0 = 0.0;
		ln_gamma_complex(1.5 - 0.5*nu_real, -0.5*nu_imag, &a_0, &b_0);
		double lgamma_num_real = 0.0; double lgamma_num_imag = 0.0;
		ln_gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);

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
		ln_gamma_complex(l + 0.5*nu_real, 0.5*nu_imag, &lgamma_num_real, &lgamma_num_imag);
		double scnd_den_gamma = ln_gamma_real(l + 1.5);
		double exp_factor = exp(lgamma_num_real - scnd_den_gamma);
		double exp_real = exp_factor*cos(lgamma_num_imag);
		double exp_imag = exp_factor*sin(lgamma_num_imag);

		//Calculate fraction of numerator and denominator
    pref_real = (first_den_gamma_real*exp_real+first_den_gamma_imag*exp_imag)/ first_den_gamma_abs_squared;
		pref_imag = (first_den_gamma_real*exp_imag-first_den_gamma_imag*exp_real)/ first_den_gamma_abs_squared;
	}
  //We have found the prefactor in an overflow-safe way
  //Now comes another factor of 2^nu
  double exp_factor = exp(LOG_2*nu_real);
  double tot_pref_real = MATH_PI*MATH_PI*exp_factor*(pref_real*cos(LOG_2*nu_imag)-pref_imag*sin(LOG_2*nu_imag));
  double tot_pref_imag = MATH_PI*MATH_PI*exp_factor*(pref_real*sin(LOG_2*nu_imag)+pref_imag*cos(LOG_2*nu_imag));
  //Now comes the actual taylor expansion
  
  //Setup summation variables
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
  double gam = 1.;
  double pow2 = 2.;
  //Do a Taylor summation up to N_BESSEL_LTAYLOR_COEFFS
  for(i=0;i<N_BESSEL_LTAYLOR_COEFFS;++i){
    taylor_coefficient = 1./(gam*pow2);//equivalent to 1.0/(gamma_real(i+1.0)*pow(2.0,i+1.0));
    gam*=(i+1.0);
    pow2*=2;
    cur_val_real = taylor_coefficient*prod_fac_real*cur_t;
    cur_val_imag = taylor_coefficient*prod_fac_imag*cur_t;
    cur_sum_real += cur_val_real;
    cur_sum_imag += cur_val_imag;
    cur_t*=t*t;

    cur_fac_real = ((-1.0+nu_real+2.0*i)*(l+0.5*nu_real+i)-0.5*nu_imag*nu_imag)/(l+1.5+i);
    cur_fac_imag = ((-1.0+nu_real+2.0*i)*0.5*nu_imag+(l+0.5*nu_real+i)*nu_imag)/(l+1.5+i);

    temp_real = (cur_fac_real*prod_fac_real-cur_fac_imag*prod_fac_imag);
    temp_imag = (cur_fac_imag*prod_fac_real+cur_fac_real*prod_fac_imag);
    prod_fac_real = temp_real;
    prod_fac_imag = temp_imag;

  }
  //Get result
  *res_real = pref_t*(tot_pref_real*cur_sum_real-tot_pref_imag*cur_sum_imag);
  *res_imag = pref_t*(tot_pref_imag*cur_sum_real+tot_pref_real*cur_sum_imag);
  return;
}
/**
 * This function calculates the integral over the bessel functions
 *
 * When l=0 or l=1 the function is a very simple analytical one
 *
 * When getting into regimes of t<<1 or (1-t*t)<<1,
 *  we have to switch the method of calculation
 * This is done using a BESSEL_INTEGRAL_T_SWITCH
 *  which decides how high t*t should be for a switch to occur
 * Additionally, BESSEL_INTEGRAL_T_SWITCH_SAFE signals,
 *  for which t the safe version can be used.
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
		bessel_integral_lowt_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
  else if(t*t<BESSEL_INTEGRAL_T_SWITCH){
    bessel_integral_lowt_transform_safe(l, nu_real, nu_imag, t, res_real,res_imag);
    return;
  }
	else{
		bessel_integral_hight_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
}
/**
 * This function calculates the integral over the bessel functions,
 *  in an overflow-safe version
 *
 * When l=0 or l=1 the function is a very simple analytical one
 *
 * When getting into regimes of t<<1 or (1-t*t)<<1,
 *  we have to switch the method of calculation
 * This is done using a BESSEL_INTEGRAL_T_SWITCH
 *  which decides how high t*t should be for a switch to occur
 * Additionally, BESSEL_INTEGRAL_T_SWITCH_SAFE signals,
 *  for which t the safe version can be used.
 * 
 * Additionally, if overflow_flag is still _TRUE_, then the 
 *  overflow-safe version will be used, despite being in the 
 *  likely safe regime. (For large nu, the safe regime might
 *  not be as "safe" as expected)
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
		bessel_integral_lowt_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
  else if(t*t<BESSEL_INTEGRAL_T_SWITCH){
    bessel_integral_lowt_transform_safe_flag(l, nu_real, nu_imag, t, res_real,res_imag,overflow_flag);
    return;
  }
	else{
		bessel_integral_hight_transform(l, nu_real, nu_imag, t, res_real,res_imag);
		return;
	}
}

/**
 * The Hypergeometric function for l = 0 is known analytically. The below implements the corresponding formula.
 *  There is a small expansion required for t<<1 due to numerical divergences of the full formula.
 * 
 * This is the HYPERGEOMETRIC function at l=0, not the BESSEL INTEGRAL
 * This is only the small helper-function required for complex backward
 *  "helper" function recursion
 * 
 * ! Never insert t=1
 * */
#define NUM_HYPERGEOM_TAYLOR_COEFF 15
#define BESSEL_HYPERGEOM_TAYLOR 0.01
void bessel_hypergeom_l0(double nu_real,double nu_imag,double t,double* res_real,double* res_imag){
  //If t<BESSEL_HYPERGEOM_TAYLOR use Taylor, otherwise full
  if(t>BESSEL_HYPERGEOM_TAYLOR){
    double log1pt = log(1.0+t);
    double log1mt = log(1.0-t); //Can be -inf for t=1.0
    //Use full formula
    double exp_factor = exp(log1pt*(2-nu_real));
    double first_real = -exp_factor*cos(-log1pt*nu_imag);
    double first_imag = -exp_factor*sin(-log1pt*nu_imag);
    exp_factor = exp(log1mt*(2-nu_real));
    double second_real = exp_factor*cos(-log1mt*nu_imag);
    double second_imag = exp_factor*sin(-log1mt*nu_imag);
    double den_abs_sq = (nu_real-2.0)*(nu_real-2.0)+nu_imag*nu_imag;
    // Get result
    *res_real = 0.5/t*((first_real+second_real)*(nu_real-2.0)+(first_imag+second_imag)*nu_imag)/den_abs_sq;
    *res_imag = 0.5/t*((first_imag+second_imag)*(nu_real-2.0)-(first_real+second_real)*nu_imag)/den_abs_sq;
    return;
  }
  else{
    //If t is too small the above depends on precise cancellations
    //Instead, we use this below, which expands in small t
    
    //Setup Talyor summation
    double cur_nu_fac_real=1.0,cur_nu_fac_imag=0.0;
    double cur_t_fac = 1.0;
    double sum_real = 0.0;
    double sum_imag = 0.0;
    double temp_real,temp_imag;
    double temp_nu_real,temp_nu_imag;
    int i;
    double gam = 1.0;
    //Use Taylor expansion at low t^2 << 1
    for(i=0;i<NUM_HYPERGEOM_TAYLOR_COEFF;++i){
      sum_real += cur_nu_fac_real*cur_t_fac/gam;
      sum_imag += cur_nu_fac_imag*cur_t_fac/gam;
      gam*=(2*i+2)*(2*i+3);
      //Calculate new term
      cur_t_fac*=t*t;
      temp_nu_real = (nu_real+2*i-1)*(nu_real+2*i)-nu_imag*nu_imag;
      temp_nu_imag = (nu_real+2*i-1)*nu_imag+(nu_real+2*i)*nu_imag;
      temp_real = (temp_nu_real*cur_nu_fac_real-temp_nu_imag*cur_nu_fac_imag);
      temp_imag = (temp_nu_real*cur_nu_fac_imag+temp_nu_imag*cur_nu_fac_real);
      //Add new term to sum
      cur_nu_fac_real = temp_real;
      cur_nu_fac_imag = temp_imag;
    }
    // Get result
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
}
/**
 * The Hypergeometric function for l = 1 is known analytically. The below implements the corresponding formula.
 *  There is a small expansion required for t<<1 due to numerical divergences of the full formula.
 * 
 * This is the HYPERGEOMETRIC function at l=1, not the BESSEL INTEGRAL
 * This is only the small helper-function required for complex backward
 *  "helper" function recursion
 * 
 * ! Never insert t=1
 * */
void bessel_hypergeom_l0_p1(double nu_real,double nu_imag,double t,double* res_real,double* res_imag){
  //If t<BESSEL_HYPERGEOM_TAYLOR use Taylor, otherwise full
  if(t>BESSEL_HYPERGEOM_TAYLOR){
    double log1pt = log(1.0+t);
    double log1mt = log(1.0-t); //Can be -inf for t=1.0
    //Use full formula
    double den_abs_sq = (nu_real-2.0)*(nu_real-2.0)+nu_imag*nu_imag;

    double exp_factor = exp(log1mt*(1-nu_real));
    double first_real = exp_factor*cos(-log1mt*nu_imag);
    double first_imag = exp_factor*sin(-log1mt*nu_imag);
    exp_factor = exp(log1pt*(1-nu_real));
    double second_real = exp_factor*cos(-log1pt*nu_imag);
    double second_imag = exp_factor*sin(-log1pt*nu_imag);

    //Take care of the first of two different brackets
    double first_bracket_real = (1-t)*first_real-(1+t)*second_real;
    double first_bracket_imag = (1-t)*first_imag-(1+t)*second_imag;
    double pref_first_bracket_real = 1/t*((nu_real-1.0)*(nu_real-2.0)+nu_imag*nu_imag)/den_abs_sq;
    double pref_first_bracket_imag = 1/t*(nu_imag*(nu_real-2.0)-(nu_real-1.0)*nu_imag)/den_abs_sq;

    double tot_first_bracket_real = first_bracket_real*pref_first_bracket_real-first_bracket_imag*pref_first_bracket_imag;
    double tot_first_bracket_imag = first_bracket_imag*pref_first_bracket_real+first_bracket_real*pref_first_bracket_imag;

    //Take care of the second of two different brackets
    double second_bracket_real = first_real+second_real;
    double second_bracket_imag = first_imag+second_imag;

    double tot_brackets_real = tot_first_bracket_real+second_bracket_real;
    double tot_brackets_imag = tot_first_bracket_imag+second_bracket_imag;

    den_abs_sq = nu_real*nu_real+nu_imag*nu_imag;
    // Get result
    *res_real = 0.5*(tot_brackets_real*nu_real+tot_brackets_imag*nu_imag)/den_abs_sq;
    *res_imag = 0.5*(tot_brackets_imag*nu_real-tot_brackets_real*nu_imag)/den_abs_sq;
    return;
  }
  else{
    //If t is too small the above depends on precise cancellations
    //Instead, we use this below, which expands in small t
    
    //Setup taylor summation
    double cur_nu_fac_real=1.0,cur_nu_fac_imag=0.0;
    double cur_t_fac = 1.0;
    double sum_real = 0.0;
    double sum_imag = 0.0;
    double temp_real,temp_imag;
    double temp_nu_real,temp_nu_imag;
    int i;
    double gam = 1.0;
    //Use talyor summation at t^2 << 1
    for(i=0;i<NUM_HYPERGEOM_TAYLOR_COEFF;++i){
      sum_real += cur_nu_fac_real*cur_t_fac/gam;
      sum_imag += cur_nu_fac_imag*cur_t_fac/gam;
      gam*=(2*i+2)*(2*i+3);

      //Calculate new term in sum
      cur_t_fac*=t*t;
      temp_nu_real = (nu_real+2*i-1)*(nu_real+2*i+2)-nu_imag*nu_imag;
      temp_nu_imag = (nu_real+2*i-1)*nu_imag+(nu_real+2*i+2)*nu_imag;
      temp_real = (temp_nu_real*cur_nu_fac_real-temp_nu_imag*cur_nu_fac_imag);
      temp_imag = (temp_nu_real*cur_nu_fac_imag+temp_nu_imag*cur_nu_fac_real);
      //Add term to sum
      cur_nu_fac_real = temp_real;
      cur_nu_fac_imag = temp_imag;
    }
    //Get result
    *res_real = sum_real;
    *res_imag = sum_imag;
    return;
  }
}

/**
 * This function calculates the integral over the bessel functions,
 *  in a safe recursion, using an additional "helper" function
 *  given as 2F1(a,b+1,c|z) for a given 2F1(a,b,c|z)
 * */
//Define local variables, where GIANT_VAL*GIANT_VAL has to be double representable
#define GIANT_VAL 1e120
#define NON_INVERTIBLE_LIMIT 1e-25
int bessel_integral_recursion_complicated(int l_max,
                                          int l_recursion_max,
                                          double nu_real,
                                          double nu_imag,
                                          double t,
                                          double bi_allowed_error,
                                          double* bi_real,
                                          double* bi_imag,
                                          double* max_t,
                                          double* initial_abs,
                                          ErrorMsg errmsg){
  //Local variables, 2x2 matrix m, and vector f of the two 2F1 functions
  double z = t*t;
  double a_real = 0.5*(nu_real-1.0);
  double a_imag = 0.5*nu_imag;
  double m11_real = 1.0,m11_imag=0.0;
  double m12_real = 0.0,m12_imag=0.0;
  double m21_real = 0.0,m21_imag=0.0;
  double m22_real = 1.0,m22_imag=0.0;
  double f000_real=1.0,f000_imag=0.0;
  double f010_real=1.0,f010_imag=0.0;
  //Temporary variables during matrix multiplication step
  double temp1_real,temp1_imag;
  double temp2_real,temp2_imag;
  double temp3_real,temp3_imag;
  double temp4_real,temp4_imag;
  double tempf_real,tempf_imag;
  //Setup the initial variables
  bi_real[l_max] = 1.0;
  bi_imag[l_max] = 0.0;
  //Loop over recursion l
  int index_l;
  for(index_l=l_recursion_max;index_l>0;--index_l){
    //Calculate local loop variables
    double c = 1.5+index_l-1;
    double alpha_real = a_real/c;
    double alpha_imag = a_imag/c;
    double z_alpha_real = z-alpha_real;
    double o_alpha_real = 1.+alpha_real;
    //Calculate new matrix

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

    tempf_real = -alpha_real*f000_real+alpha_imag*f000_imag+(o_alpha_real)*f010_real-alpha_imag*f010_imag;
    tempf_imag = -alpha_real*f000_imag-alpha_imag*f000_real+(o_alpha_real)*f010_imag+alpha_imag*f010_real;

    //Prepare (f000,f010) vector for next step
    f010_real = tempf_real;
    f010_imag = tempf_imag;
    f000_real = bi_real[index_l-1];
    f000_imag = bi_imag[index_l-1];

    //Prepare matrix for next step
    m11_real = temp1_real;
    m11_imag = temp1_imag;
    m12_real = temp2_real;
    m12_imag = temp2_imag;
    m21_real = temp3_real;
    m21_imag = temp3_imag;
    m22_real = temp4_real;
    m22_imag = temp4_imag;
  }

  //Now we know the backward recursion-vector starting from from (1+0i,1+0i), 
  // including all bessel values inbetween
  // and also the transition matrix.
  double ad_real = (m11_real*m22_real-m11_imag*m22_imag);
  double ad_imag = (m11_real*m22_imag+m11_imag*m22_real);
  double bc_real = (m12_real*m21_real-m12_imag*m21_imag);
  double bc_imag = (m12_real*m21_imag+m12_imag*m21_real);

  //Check for invertibility of transition matrix, by first checking for overflows during check
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

  //Now do the actual invertibility check
  double det_abs = (ad_real-bc_real)*(ad_real-bc_real)+(ad_imag-bc_imag)*(ad_imag-bc_imag);
  double ad_abs = ad_real*ad_real+ad_imag*ad_imag;
  double bc_abs = bc_real*bc_real+bc_imag*bc_imag;

  //If this matrix is invertible, we are in the forward recursion case
  //Otherwise, we have already found the results by multiplying with lambda
  
  //Compare to analytical result at l=0
  double res_real,res_imag;
  bessel_hypergeom_l0(nu_real,nu_imag,t,&res_real,&res_imag);
  //If invertible...
  if(det_abs/(ad_abs+bc_abs) < NON_INVERTIBLE_LIMIT){
    //Do the inversion to find lambda, and multiply all factors
    double lambda_real = (f000_real*res_real+f000_imag*res_imag)/(res_real*res_real+res_imag*res_imag);
    double lambda_imag = (f000_imag*res_real-f000_real*res_imag)/(res_real*res_real+res_imag*res_imag);

    double lambda_abs_sq = lambda_real*lambda_real+lambda_imag*lambda_imag;

    double temp_real,temp_imag;
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      temp_real = (bi_real[index_l]*lambda_real+bi_imag[index_l]*lambda_imag)/lambda_abs_sq;
      temp_imag = (bi_imag[index_l]*lambda_real-bi_real[index_l]*lambda_imag)/lambda_abs_sq;
      bi_real[index_l] = temp_real;
      bi_imag[index_l] = temp_imag;
    }
    //Check for the error. Impressively, this is actually possible here.
    double err = sqrt(((bi_real[0]-res_real)*(bi_real[0]-res_real)+(bi_imag[0]-res_imag)*(bi_imag[0]-res_imag))/(res_real*res_real+res_imag*res_imag));
    if(err>bi_allowed_error){
      sprintf(errmsg,"%s(L:%d) : Backwards recursion for bessel integrals returned unnaturally high error. Outside of acceptable parameter region assumed.",__func__,__LINE__);
      return _FAILURE_;
    }
  }
  //If not invertible...
  else{
    //Forward recursion should be stable... Set up forward initial conditions
    bi_real[0] = res_real;
    bi_imag[0] = res_imag;
    double temp1_real,temp1_imag;
    double temp2_real,temp2_imag;
    double begin_real,begin_imag;
    bessel_hypergeom_l0_p1(nu_real,nu_imag,t,&begin_real,&begin_imag);
    //Start forward recursion
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      //Local loop variables
      double c = 1.5+index_l;
      double alpha_real = a_real/c;
      double alpha_imag = a_imag/c;

      double alpha_m1_abs_sq = (1-alpha_real)*(1-alpha_real)+alpha_imag*alpha_imag;
      double alpha_m1_sq_real = 1.0-alpha_real*alpha_real+alpha_imag*alpha_imag;
      double alpha_m1_sq_imag =    -alpha_real*alpha_imag-alpha_imag*alpha_real;

      //Calculate transition matrix in forward direction
      double den_abs_sq = alpha_m1_sq_real*alpha_m1_sq_real+alpha_m1_sq_imag*alpha_m1_sq_imag;
      double c_real = (alpha_real*alpha_m1_sq_real+alpha_imag*alpha_m1_sq_imag)/den_abs_sq/z;
      double c_imag = (alpha_imag*alpha_m1_sq_real-alpha_real*alpha_m1_sq_imag)/den_abs_sq/z;
      double d_real = ( (z-alpha_real)*alpha_m1_sq_real-alpha_imag*alpha_m1_sq_imag)/den_abs_sq/z;
      double d_imag = (-(z-alpha_real)*alpha_m1_sq_imag-alpha_imag*alpha_m1_sq_real)/den_abs_sq/z;

      //Calculate next terms in forward direction
      temp1_real = ((1-alpha_real)*res_real-alpha_imag*res_imag-(1-z)*(1-alpha_real)*begin_real+(1-z)*alpha_imag*begin_imag)/alpha_m1_abs_sq/z;
      temp1_imag = ((1-alpha_real)*res_imag+alpha_imag*res_real-(1-z)*(1-alpha_real)*begin_imag-(1-z)*alpha_imag*begin_real)/alpha_m1_abs_sq/z;
      temp2_real = (c_real*res_real-c_imag*res_imag+d_real*begin_real-d_imag*begin_imag);
      temp2_imag = (c_real*res_imag+c_imag*res_real+d_real*begin_imag+d_imag*begin_real);

      //Store values
      bi_real[index_l+1] = temp1_real;
      bi_imag[index_l+1] = temp1_imag;
      begin_real = temp2_real;
      begin_imag = temp2_imag;
      res_real = bi_real[index_l+1];
      res_imag = bi_imag[index_l+1];
    }
  }
  //The bi_real and bi_imag arrays are now filled with Hypergeometric2F1((nu-1)/2,l+nu/2,l+3/2,t*t),
  // but now we want to include the prefactors
  double temp_real,temp_imag;
  double exp_factor = MATH_PI*MATH_PI*exp(LOG_2*(nu_real-1.0));
  double pref_pref_real = exp_factor*cos(LOG_2*nu_imag);
  double pref_pref_imag = exp_factor*sin(LOG_2*nu_imag);
  double start_frac_real,start_frac_imag;
  //Overflow safe version of getting the INITIAL Gamma functions for l=0
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
    ln_gamma_complex(1.5-0.5*nu_real,-0.5*nu_imag,&a0,&b0);
    ln_gamma_complex(0.5*nu_real,0.5*nu_imag,&a1,&b1);
    exp_factor = exp(a1-a0);
    start_frac_real = exp_factor*cos(b1-b0);
    start_frac_imag = exp_factor*sin(b1-b0);
  }
  temp_real = 2.0/sqrt(MATH_PI)*(pref_pref_real*start_frac_real-pref_pref_imag*start_frac_imag);
  temp_imag = 2.0/sqrt(MATH_PI)*(pref_pref_imag*start_frac_real+pref_pref_real*start_frac_imag);
  start_frac_real= temp_real;
  start_frac_imag= temp_imag;

  //Obtain all further Gamma functions by their recursion relation Gamma(1+z)=z*Gamma(z)
  double cur_frac_real = start_frac_real;
  double cur_frac_imag = start_frac_imag;
  for(index_l=0;index_l<=l_max;++index_l){
    temp_real = cur_frac_real*bi_real[index_l]-cur_frac_imag*bi_imag[index_l];
    temp_imag = cur_frac_imag*bi_real[index_l]+cur_frac_real*bi_imag[index_l];
    bi_real[index_l] = temp_real;
    bi_imag[index_l] = temp_imag;

    //If the integral has become very small, allow it to "exit" the recursion relation
    if(bi_real[index_l]*bi_real[index_l]+bi_imag[index_l]*bi_imag[index_l]<BESSEL_EPSILON*initial_abs[index_l]){
      max_t[index_l] = t;
    }
    temp_real = t*(cur_frac_real*(index_l+0.5*nu_real)-cur_frac_imag*0.5*nu_imag)/(index_l+1.5);
    temp_imag = t*(cur_frac_imag*(index_l+0.5*nu_real)+cur_frac_real*0.5*nu_imag)/(index_l+1.5);
    cur_frac_real = temp_real;
    cur_frac_imag = temp_imag;
  }
  return _SUCCESS_;
}

/**
 * This function calculates the integral over the bessel functions
 *  at t=1 by the simple forward recursion relation
 * 
 * Explicitly t=1 is important for criterium of 'vanishing' function
 *  (as a good approximation for maximum value for t in [0,1])
 * */
int bessel_integral_recursion_initial_abs(int l_max,double nu_real,double nu_imag,double* abi_real,double* abi_imag,double* initial_abs){
  //Initialize recursion relation
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

  //Initialize loop
  double bi_next_next_real;
  double bi_next_next_imag;
  double temp_real, temp_imag;
  int index_l;
  for(index_l=0;index_l<=l_max;++index_l){
    //Store loop variables
    double l = (double)index_l;
    abi_real[index_l]= bi_real;
    abi_imag[index_l]= bi_imag;
    initial_abs[index_l]=bi_real*bi_real+bi_imag*bi_imag;
    //Calculate next loop iteration using forward recursion
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
    //Add term to recursion relation
    bi_real = bi_next_real;
    bi_imag = bi_next_imag;
    bi_next_real = bi_next_next_real;
    bi_next_imag = bi_next_next_imag;
  }
  return _SUCCESS_;
}

/**
 * This function calculates the integral over the bessel functions
 *  at small t by a simple recursion relation
 * */
int bessel_integral_recursion_taylor(int l_max,double nu_real,double nu_imag,double t,double* max_t,double* initial_abs,double* bi_real,double* bi_imag){
  double res_real,res_imag;
  int index_l;
  for(index_l=0;index_l<=l_max;++index_l){
    //If the mode has already exited, ignore it
    if(t<max_t[index_l]){
      bi_real[index_l]=0.0;
      bi_imag[index_l]=0.0;
    }
    //Otherwise do full calculation using the talyor already implemented
    else{
      double l = (double)index_l;
      bessel_integral_transform(l,nu_real,nu_imag,t,&res_real,&res_imag);
      //Finally check for mode exiting
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
  return _SUCCESS_;
}


/**
 * This function calculates the integral over the bessel functions
 *  by recursion in the backwards direction using an overflow-safe algorithm
 * */
int bessel_integral_recursion_backward_simple_safe(
                                             int l_max,
                                             int l_recursion_max,
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             short* overflow_flag,
                                             ErrorMsg errmsg
                                            ){
  // Setup initial variables
  double bi_real,bi_imag,bi_next_real,bi_next_imag,bi_next_next_real,bi_next_next_imag;
  double temp_real,temp_imag;
  int l_start =l_recursion_max;
  //Get seeding values of Il
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

  //Store seeding values and prepare loop
  abi_real[l_start] = bi_real;
  abi_imag[l_start] = bi_imag;
  abi_real[l_start-1] = bi_next_real;
  abi_imag[l_start-1] = bi_next_imag;

  //Loop in backward direction
  int index_l;
  for(index_l=l_start;index_l>=2;--index_l){
    double l = (double)index_l;

    //Calculate with simple recursion relation the transition factors
    double den_real = (l-2.0+nu_real*0.5);
    double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
    double bi_factor_real = -((1.0+l-nu_real*0.5)*bi_real+nu_imag*bi_imag*0.5);
    double bi_factor_imag = -((1.0+l-nu_real*0.5)*bi_imag-nu_imag*bi_real*0.5);
    double bi_factor_next_real = (1+t*t)/t*(l-0.5)*bi_next_real;
    double bi_factor_next_imag = (1+t*t)/t*(l-0.5)*bi_next_imag;

    //Transition all bessel integrals
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
  }
  //End l


  //TODO :: evaluate, if there is a simple criterion when to do 
  // back-checking or not. For now, we always do it.
  int DOES_BACK_CHECK = _TRUE_;
  if(DOES_BACK_CHECK){
    //Get value at l=0
    double res_real,res_imag;
    bessel_integral_transform(0,nu_real,nu_imag,t,&res_real,&res_imag);
    double abi_0_abs_sq = abi_real[0]*abi_real[0]+abi_imag[0]*abi_imag[0];

    //Get multiplicative factor
    double lambda_real = (res_real*abi_real[0]+res_imag*abi_imag[0])/abi_0_abs_sq;
    double lambda_imag = (res_imag*abi_real[0]-res_real*abi_imag[0])/abi_0_abs_sq;

    //And finally multiply all values by multiplicative factor
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      temp_real = abi_real[index_l]*lambda_real-abi_imag[index_l]*lambda_imag;
      temp_imag = abi_real[index_l]*lambda_imag+abi_imag[index_l]*lambda_real;
      abi_real[index_l] = temp_real;
      abi_imag[index_l] = temp_imag;
      //Take care of modes exiting the recursion relation
      if(temp_real*temp_real+temp_imag*temp_imag<BESSEL_EPSILON*initial_abs[index_l]){
        max_t[index_l]=t;
      }
    }
  }
  return _SUCCESS_;
}



/**
 * This function calculates the integral over the bessel functions
 *  by recursion in the backwards direction
 * 
 * ! NOT overflow-safe
 * */
int bessel_integral_recursion_backward_simple(
                                             int l_max,
                                             int l_recursion_max,
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             ErrorMsg errmsg
                                            ){
  // Setup initial variables
  double bi_real,bi_imag,bi_next_real,bi_next_imag,bi_next_next_real,bi_next_next_imag;
  double temp_real,temp_imag;
  //Get seeding values of Il
  int l_start = l_recursion_max;
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

  //Store seeding values and prepare loop
  abi_real[l_start] = bi_real;
  abi_imag[l_start] = bi_imag;
  abi_real[l_start-1] = bi_next_real;
  abi_imag[l_start-1] = bi_next_imag;

  //Loop in backward direction
  int index_l;
  for(index_l=l_start;index_l>=2;--index_l){
    double l = (double)index_l;

    //Calculate with simple recursion relation the transition factors
    double den_real = (l-2.0+nu_real*0.5);
    double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
    double bi_factor_real = -((1.0+l-nu_real*0.5)*bi_real+nu_imag*bi_imag*0.5);
    double bi_factor_imag = -((1.0+l-nu_real*0.5)*bi_imag-nu_imag*bi_real*0.5);
    double bi_factor_next_real = (1+t*t)/t*(l-0.5)*bi_next_real;
    double bi_factor_next_imag = (1+t*t)/t*(l-0.5)*bi_next_imag;

    //Transition all bessel integrals
    temp_real = (bi_factor_next_real+bi_factor_real);
    temp_imag = (bi_factor_next_imag+bi_factor_imag);
    bi_next_next_real = (temp_real*den_real+temp_imag*nu_imag*0.5)/den_abs_squared;
    bi_next_next_imag = (temp_imag*den_real-temp_real*nu_imag*0.5)/den_abs_squared;
    bi_real = bi_next_real;
    bi_imag = bi_next_imag;
    bi_next_real = bi_next_next_real;
    bi_next_imag = bi_next_next_imag;

    //Store new integrals
    abi_real[index_l-2] = bi_next_next_real;
    abi_imag[index_l-2] = bi_next_next_imag;
  }
  //End l


  //TODO :: evaluate, if there is a simple criterion when to do 
  // back-checking or not. For now, we always do it.
  int DOES_BACK_CHECK = _TRUE_;
  if(DOES_BACK_CHECK){
    //Get value at l=0
    double res_real,res_imag;
    bessel_integral_transform(0,nu_real,nu_imag,t,&res_real,&res_imag);
    double abi_0_abs_sq = abi_real[0]*abi_real[0]+abi_imag[0]*abi_imag[0];

    //Get multiplicative factor
    double lambda_real = (res_real*abi_real[0]+res_imag*abi_imag[0])/abi_0_abs_sq;
    double lambda_imag = (res_imag*abi_real[0]-res_real*abi_imag[0])/abi_0_abs_sq;
    

    //And finally multiply all values by multiplicative factor
    int index_l;
    for(index_l=0;index_l<=l_max;++index_l){
      temp_real = abi_real[index_l]*lambda_real-abi_imag[index_l]*lambda_imag;
      temp_imag = abi_real[index_l]*lambda_imag+abi_imag[index_l]*lambda_real;
      abi_real[index_l] = temp_real;
      abi_imag[index_l] = temp_imag;
      //Take care of modes exiting the recursion relation
      if(temp_real*temp_real+temp_imag*temp_imag<BESSEL_EPSILON*initial_abs[index_l]){
        max_t[index_l]=t;
      }
    }
  }
  return _SUCCESS_;
}


/**
 * This function calculates the integral over the bessel functions
 *  by recursion in the backwards direction
 * 
 * ! NOT overflow-safe in theory, but practically not yet ever occured
 * */
int bessel_integral_recursion_forward_simple(int l_max,
                                             double nu_real,
                                             double nu_imag,
                                             double t,
                                             double* abi_real,
                                             double* abi_imag,
                                             double* max_t,
                                             double* initial_abs,
                                             ErrorMsg errmsg
                                             ){
  // Setup initial variables
  double temp_real,temp_imag;
  double bi_real,bi_imag,bi_next_real,bi_next_imag,bi_next_next_real,bi_next_next_imag;
  int l_start = 0;

  // Get seeding values
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

  //Store initial values
  abi_real[l_start] = bi_real;
  abi_imag[l_start] = bi_imag;
  abi_real[l_start+1] = bi_next_real;
  abi_imag[l_start+1] = bi_next_imag;

  //Loop in forward direction
  int index_l;
  for(index_l=l_start;index_l<=l_max;++index_l){
      //Store new variables
      abi_real[index_l] = bi_real;
      abi_imag[index_l] = bi_imag;

      //Take care of exiting modes
      if(bi_real*bi_real+bi_imag*bi_imag<BESSEL_EPSILON*initial_abs[index_l]){
        max_t[index_l]=t;
      }
      double l = (double)index_l;

      //Calculate transition factors
      double den_real = (3.0+l-nu_real*0.5);
      double den_abs_squared = den_real*den_real+0.25*nu_imag*nu_imag;
      double bi_factor_real = -((l+nu_real*0.5)*bi_real-nu_imag*bi_imag*0.5);
      double bi_factor_imag = -((l+nu_real*0.5)*bi_imag+nu_imag*bi_real*0.5);
      double bi_factor_next_real = (1+t*t)/t*(l+1.5)*bi_next_real;
      double bi_factor_next_imag = (1+t*t)/t*(l+1.5)*bi_next_imag;
      temp_real = (bi_factor_next_real+bi_factor_real);
      temp_imag = (bi_factor_next_imag+bi_factor_imag);

      //Advance Bessel integrals
      bi_next_next_real = (temp_real*den_real-temp_imag*nu_imag*0.5)/den_abs_squared;
      bi_next_next_imag = (temp_imag*den_real+temp_real*nu_imag*0.5)/den_abs_squared;
      bi_real = bi_next_real;
      bi_imag = bi_next_imag;
      bi_next_real = bi_next_next_real;
      bi_next_imag = bi_next_next_imag;
  }
  //End l
  return _SUCCESS_;
}

/**
 * This function calculates the integral over the bessel functions
 *  by tayloring in (1-t*t)
 * 
 * */
int bessel_integral_recursion_inverse_self(int l_max,
                                           double nu_real,
                                           double nu_imag,
                                           double t,
                                           double* abi_real,
                                           double* abi_imag,
                                           double* max_t,
                                           double* initial_abs,
                                           ErrorMsg errmsg){
  //Setup local variables
  double res_real,res_imag;
  double pref_first_real,pref_first_imag;
  double pref_scnd_real,pref_scnd_imag;
  double tot_pref_real,tot_pref_imag;
  double temp_real,temp_imag;
  double div_real,div_imag;
  //Common factors
  double z = t*t;
  double z_fac1 = ((1-z)*(1-z))/((1+z)*(1+z));
  double z_fac2 = 2*(1+z)/(1-z);
  double lnz_fac3 = log(2.0/(1.0+z));
  double l = 0;
  //Overflow safe-versions for INITIAL prefactors at l=0
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
    ln_gamma_complex(1.5-nu_real*0.5, -nu_imag*0.5, &a_0, &b_0);
    double a_1 = 0.0; double b_1 = 0.0;
    ln_gamma_complex(l + nu_real*0.5, nu_imag*0.5, &a_1, &b_1);
    double a_2 = 0.0; double b_2 = 0.0;
    ln_gamma_complex(1 - nu_real*0.5, -nu_imag*0.5, &a_2, &b_2);
    double a_3 = 0.0; double b_3 = 0.0;
    ln_gamma_complex(l + 2 - nu_real*0.5, -nu_imag*0.5, &a_3, &b_3);

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
    ln_gamma_complex(nu_real*0.5 - 1, nu_imag*0.5, &c_0, &d_0);
    exp_factor = exp(c_0-a_0);
    exp_real = exp_factor*cos(d_0-b_0);
    exp_imag = exp_factor*sin(d_0-b_0);

    //Calculate second prefactor
    double pref_num_gamma_real = exp_real;
    double pref_num_gamma_imag = exp_imag;
    pref_scnd_real = pre_pref_scnd_real*pref_num_gamma_real - pre_pref_scnd_imag*pref_num_gamma_imag;
    pref_scnd_imag = pre_pref_scnd_real*pref_num_gamma_imag + pre_pref_scnd_imag*pref_num_gamma_real;
  }


  //Now loop through all l and each time get the Taylor approximation
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

    //Update factors
    tot_pref_real*=(2./(1.+z))*t;
    tot_pref_imag*=(2./(1.+z))*t;

    div_real = l+1-nu_real/2. +1;
    div_imag = -nu_imag/2.;
    double div_abs_sq = div_real*div_real+div_imag*div_imag;
    double factor_real = ((l-1+nu_real*0.5 +1)*div_real+div_imag*nu_imag*0.5)/div_abs_sq;
    double factor_imag = (nu_imag*0.5*div_real-(1+   l-1+nu_real*0.5)*div_imag)/div_abs_sq;

    //Advance values
    temp_real = factor_real*pref_first_real-factor_imag*pref_first_imag;
    temp_imag = factor_real*pref_first_imag+factor_imag*pref_first_real;

    pref_first_real = temp_real;
    pref_first_imag = temp_imag;

    //Check for modes "exiting" the recursion relations
    if(res_real*res_real+res_imag*res_imag<BESSEL_EPSILON*initial_abs[index_l]){
      max_t[index_l]=t;
    }

    //Store values of the Bessel integrals
    abi_real[index_l] = res_real;
    abi_imag[index_l] = res_imag;
  }
  return _SUCCESS_;
}
