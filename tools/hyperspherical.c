/** @file hyperspherical.c Documented hyperspherical bessel function module.
 *
 * Thomas Tram, 11.01.2013
 *
 * This module computes hyperspherical Bessel functions.
 */

#include "hyperspherical.h"

int HypersphericalClosedGegenbauer(int l, int beta, double y, double *Phi, ErrorMsg errmsg){
/** Calculate hyperspherical function using Gegenbauer polynomials
 */
  double MM;
  //Calculate MM avoiding overflow:
  if (l<80){
    MM = pow(2,l)*Hypfactorial(l)*
      sqrt(Hypfactorial(beta-l-1)/(beta*Hypfactorial(beta+l)));
  }
  else{
    MM = exp(l*log(2)+lgamma(l+1)+0.5*lgamma(beta-l)
	     -0.5*log(beta)-0.5*lgamma(beta+l+1));
  }
  *Phi = MM*pow(sin(y),l)*Gegenbauer(beta-l-1,l+1,cos(y));

  return _SUCCESS_;
}

int HypersphericalOpenRecurrence(int l, double beta, double x, double *Phi_l, ErrorMsg errmsg){
  double maxiter_for_CF;
  double xtp_tilde;
  
  maxiter_for_CF = 100.0+2*l;
  xtp_tilde = asinh((l+maxiter_for_CF)/beta);

  if (x<xtp_tilde){
    //Use backwards recurrence from continued fraction
    HypersphericalBackward(-1, l, beta, x, Phi_l);
  }
  else{
    //Use forward recurrence:
    HypersphericalForward(-1, l, beta, x, Phi_l);
  }

  return _SUCCESS_;
}

int HypersphericalWKB(int K,int l,double beta,double y, double *Phi, ErrorMsg errmsg){
  double e, w, w2, alpha, alpha2, CscK=0, ytp, t;
  double S=0, Q, C, argu, Ai;
  int airy_sign = 1, phisign = 1, intbeta;
  double ldbl = l;

  if (K==1){
    //Limit range to [0; pi/2]:
    intbeta = (int)(beta+0.5); //Round to nearest integer (just to be sure)
    phisign =  ClosedModY(l, intbeta, &y);
  }
  e = 1.0/sqrt(ldbl*(ldbl+1.0));
  alpha = beta*e;
  if (K==-1){
    CscK = 1.0/sinh(y);
    ytp = asinh(1.0/alpha);
  }
  else if (K==1){
    CscK = 1.0/sin(y);
    ytp = asin(1.0/alpha);
  }
  else{
    //Failure
  }
  w = alpha/CscK;
  w2 = w*w;
  alpha2 = alpha*alpha;
  if (K==-1){
    if (y>ytp){
      S = alpha*log((sqrt(w2-1.0)+sqrt(w2+alpha2))/sqrt(1.0+alpha2))+
	atan(1.0/alpha*sqrt((w2+alpha2)/(w2-1.0)))-0.5*_PI_;
      airy_sign = -1;
    }
    else{
      t = sqrt(1.0-w2)/sqrt(1.0+w2/alpha2);
      S = atanh(t)-alpha*atan(t/alpha);
      airy_sign = 1;
    }
  }
  else if (K==1){
    if (y>ytp){
      t = sqrt(1-w2/alpha2)/sqrt(w2-1.0);
        S = atan(t)+alpha*atan(1.0/(t*alpha))-0.5*_PI_;
        airy_sign = -1;
    }
    else{
      S = atanh(sqrt(1.0-w2)/sqrt(1.0-w2/alpha2))-
	alpha*log((sqrt(alpha2-w2)+sqrt(1.0-w2))/sqrt(alpha2-1.0));
      airy_sign = 1;
    }
  }
  argu = 3.0*S/(2.0*e);
  Q = CscK*CscK-alpha2;
  C = 0.5*sqrt(alpha)/beta;
  Ai = airy_cheb_approx(airy_sign*pow(argu,2.0/3.0));
  *Phi = phisign*2.0*sqrt(_PI_)*C*pow(argu,1.0/6.0)*pow(fabs(Q),-0.25)*Ai*CscK;

  return _SUCCESS_;
}

int HypersphericalExplicit(int K,int l, double beta,double x, double *Phi, ErrorMsg errmsg){
  /** Explicit formulae for the Hyperspherical Besselfunctions of order l<=9.
      phi_tilde = gam * beta * cos(x*beta) + delta * sin(x*beta),
      and Phi = phi_tilde *cscK/sqrt(NK). Gamma and delta are polynomials in
      beta and cscK, containing only even powers.
  */
  double NK,xbeta,gamma,delta,CotK;
  double beta2,beta4,beta6,beta8,beta12,beta16;
  double CscK,CscK2,CscK4,CscK6,CscK8;
  //NK = prod(beta^2-K*(0:l).^2);
  beta2 = beta*beta;
  xbeta = x*beta;
  if (K==-1){
    CotK = 1.0/tanh(x);
    CscK = 1.0/sinh(x);
  }
  else if (K==1){
    CotK = 1.0/tan(x);
    CscK = 1.0/sin(x);
  }
  else{
    CotK = 1.0/x;
    CscK = CotK;
  }

  //Calculate polynomials:
  switch (l){
  case 0:
    gamma = 0;
    delta = 1;
    NK = beta2;
    break;
  case 1:
    gamma = -1;
    delta = CotK;
    NK = beta2*(beta2 - K);
    break;
  case 2:
    beta4 = beta2*beta2;
    CscK2 =CscK*CscK;
    gamma = -3*CotK;
    delta = -beta2 + 3*CscK2 - 2*K;
    NK = beta2*(4.0 + beta4 - 5.0*beta2*K);
    break;
  case 3:
    beta4 = beta2*beta2;
    CscK2 =CscK*CscK;
    gamma = beta2-15*CscK2+11*K;
    delta = CotK*(-6*beta2 + 15*CscK2 - 6*K);
    NK = beta2*(49*beta2 + beta2*beta4 - 36*K - 14*beta4*K);
    break;
  case 4:
    beta4 = beta2*beta2;
    CscK2 = CscK*CscK; CscK4 = CscK2*CscK2;
    gamma = CotK*(10*beta2-105*CscK2+50*K);
    delta = 24 + beta4 + 105*CscK4 + CscK2*(-45*beta2 - 120*K) + 35*beta2*K;
    NK = beta2*(576 + 273*beta4 + beta4*beta4 - 10*beta2*(82 + 3*beta4)*K);
    break;
  case 5:
    beta2 = beta*beta; beta4 = beta2*beta2;
    CscK2 = CscK*CscK; CscK4 = CscK2*CscK2;
    gamma = -274-beta4+105*beta2*CscK2-945*CscK4-85*beta2*K+1155*CscK2*K;
    delta = CotK*(120 + 15*beta4 + 945*CscK4 + CscK2*(-420*beta2 - 840*K) + 
		  225*beta2*K);
    NK = beta2*(beta2*(21076.0 + 1023*beta4 + beta4*beta4) - 
		5.0*(2880.0 + 11*beta4*(139 + beta4))*K);
    break;
  case 6:
    beta2 = beta*beta; beta4 = beta2*beta2; beta6 = beta4*beta2; beta8 = beta4*beta4;
    CscK2 = CscK*CscK; CscK4 = CscK2*CscK2; CscK6 = CscK4*CscK2;
    gamma = CotK*(-1764 - 21*beta4 + 1260*beta2*CscK2 - 10395*CscK4 - 
		  735*beta2*K + 10080*CscK2*K);
    delta = -1624*beta2 - beta6 + 10395*CscK6 + CscK4*(-4725*beta2 - 17010*K) - 
      720*K - 175*beta4*K + CscK2*(7560 + 210*beta4 + 6090*beta2*K);
    NK = beta2*(518400 + beta8*beta4 + 296296*beta4 + 3003*beta8 - 
		13*beta2*(59472 + 3421*beta4 + 7*beta8)*K);
    break;  
  case 7:
    beta2 = beta*beta; beta4 = beta2*beta2; beta6 = beta4*beta2; beta8 = beta4*beta4;
    beta12 = beta6*beta6;
    CscK2 = CscK*CscK; CscK4 = CscK2*CscK2; CscK6 = CscK4*CscK2;
    gamma = 6769*beta2 + beta6 - 112392*CscK2 - 378*beta4*CscK2 + 
      17325*beta2*CscK4 - 135135*CscK6 + 13068*K + 322*beta4*K - 
      23310*beta2*CscK2*K + 232155*CscK4*K;
    delta = CotK*(-13132*beta2 - 28*beta6 + 135135*CscK6 + 
		  CscK4*(-62370*beta2 - 187110*K) - 5040*K - 1960*beta4*K + 
		  CscK2*(68040 + 3150*beta4 + 64890*beta2*K));
    NK = beta2*(beta2*(38402064 + beta6*beta6 + 2475473*beta4 + 7462*beta8) - 
		20*(1270080 + 7*beta12 + 764582*beta4 + 9581*beta8)*K);
    break;  
  case 8:
    beta2 = beta*beta; beta4 = beta2*beta2; beta6 = beta4*beta2; beta8 = beta4*beta4;
    beta12 = beta6*beta6; beta16 = beta8*beta8;
    CscK2 = CscK*CscK; CscK4 = CscK2*CscK2; CscK6 = CscK4*CscK2; CscK8 = CscK4*CscK4;
    gamma = CotK*(67284*beta2 + 36*beta6 - 1191960*CscK2 - 6930*beta4*CscK2 + 
		  270270*beta2*CscK4 - 2027025*CscK6 + 109584*K + 4536*beta4*K - 
		  297990*beta2*CscK2*K + 2972970*CscK4*K);
    delta = 40320 + 22449*beta4 + beta8 + 2027025*CscK8 + 
      CscK6*(-945945*beta2 - 4324320*K) + 118124*beta2*K + 546*beta6*K + 
      CscK4*(2993760 + 51975*beta4 + 1694385*beta2*K) + 
      CscK2*(-879480*beta2 - 630*beta6 - 725760*K - 72450*beta4*K);
    NK = beta2*(1625702400 + 16422*beta12 + beta16 + 1017067024*beta4 + 14739153*beta8 - 
		68*beta2*(36516672 + 3*beta12 + 2554734*beta4 + 9841*beta8)*K);
  break;
  case 9:
    beta2 = beta*beta; beta4 = beta2*beta2; beta6 = beta4*beta2; beta8 = beta4*beta4;
    beta12 = beta6*beta6; beta16 = beta8*beta8;
    CscK2 = CscK*CscK; CscK4 = CscK2*CscK2; CscK6 = CscK4*CscK2; CscK8 = CscK4*CscK4;
    gamma = -1026576 - 63273*beta4 - beta8 + 4830210*beta2*CscK2 + 
      990*beta6*CscK2 - 55945890*CscK4 - 135135*beta4*CscK4 + 
      4729725*beta2*CscK6 - 34459425*CscK8 - 723680*beta2*K - 870*beta6*K + 
      14933160*CscK2*K + 194040*beta4*CscK2*K - 8783775*beta2*CscK4*K + 
      76351275*CscK6*K;
    delta = CotK*(362880 + 269325*beta4 + 45*beta8 + 34459425*CscK8 + 
		  CscK6*(-16216200*beta2 - 64864800*K) + 1172700*beta2*K + 9450*beta6*K + 
		  CscK4*(38918880 + 945945*beta4 + 24999975*beta2*K) + 
		  CscK2*(-10866240*beta2 - 13860*beta6 - 7983360*K - 1094940*beta4*K));
    NK = beta2*(beta2*(202759531776. + 32946*beta12 + beta16 + 15088541896.*beta4 + 68943381*beta8) - 
		5*(26336378880. + 19*beta4*(893321712 + 3*beta12 + 14395719*beta4 + 21046*beta8))*K);
    break;
  default:
    *Phi = 0.0;
    //Failure
    return _FAILURE_;
  }
  beta2 = beta*beta;
  NK = beta*beta;
  int n;
  for (n=1; n<=l; n++)
    NK *=(beta*beta-K*n*n);

  *Phi = (gamma*beta*cos(xbeta)+delta*sin(xbeta))*CscK/sqrt(NK);
  return _SUCCESS_;
}



/**
 * Gegenbauer polynomials used by closed model calculation.
 */
double Gegenbauer(int n,int alpha, double x){
  double G=0,Gm1,Gm2;
  int k;
  switch(n){
  case 0:
    G = 1.0;
    return G;
  case 1:
    G = 2*alpha*x;
    return G;
  case 2:
    G = -alpha + 2*alpha*(1+alpha)*x*x;
    return G;
  case 3:
    G = -2*alpha*(1+alpha)*x+4.0/3.0*alpha*(1+alpha)*(2+alpha)*x*x*x;
    return G;
  default:
        Gm2 = -alpha + 2*alpha*(1+alpha)*x*x;
        Gm1 = -2*alpha*(1+alpha)*x+4.0/3.0*alpha*(1+alpha)*(2+alpha)*x*x*x;
        for(k=4;k<=n; k++){
            G = (2.0*(k+alpha-1.0)*x*Gm1 - (k+2.0*alpha-2.0)*Gm2) / k;
            Gm2 = Gm1;
            Gm1 = G;
	}
	return G;
  }
}

double Hypfactorial(int n){
/** Factorial for HypersphericalClosedGegenbauer.
    It is already checked that no input creates overflow.
*/
  int j;
  double res = 1.0;
  for (j=2; j<=n; j++){
    res *= j;
  }
  return res;
}



int get_CF1(int K,int l,double beta,double x,double *CF,int *isign){
int maxiter = 100000;
 double tiny = 1e-100;
 double reltol = DBL_EPSILON;
 double CotK,aj,bj,fj,Cj,Dj,Delj;
 double beta2 = beta*beta;
 double sqrttmp;
 int j;

 if (K==-1){
    CotK = 1/tanh(x);
 }
 else if(K==1){
    CotK = 1/tan(x);
 }
 else{
    CotK = 1/x;
 }
 bj = l*CotK; //This is b_0
 fj = bj;
 Cj = bj;
 Dj = 0.0;
 *isign = 1;
 for(j=1; j<=maxiter; j++){
   sqrttmp = sqrt(beta2-K*(l+j+1)*(l+j+1));
   aj = -sqrt(beta2-K*(l+j)*(l+j))/sqrttmp;
   if (j==1)
     aj = sqrt(beta2-K*(l+1)*(l+1))*aj;
   bj = (2*(l+j)+1)/sqrttmp*CotK;
   Dj = bj+aj*Dj;
   if (Dj==0.0)
     Dj = tiny;
   Cj = bj+aj/Cj;
   if (Cj==0.0)
        Cj = tiny; 
   Dj = 1.0/Dj;
   Delj = Cj*Dj;
   fj = fj*Delj;
   if (Dj<0)
        *isign *= -1;
   if (fabs(Delj-1.0)<reltol){
     *CF = fj;
     return _SUCCESS_;
   }
 }
 return _FAILURE_;
}


int HypersphericalForward(int K,int l, double beta,double x, double *Phi_l){
  double beta2 = beta*beta;
  double CotK,Phim2,Phim1,Phi;
  int ll;

  if (K==-1){
    CotK = 1.0/tanh(x);
    Phim2 = 1.0/beta*sin(beta*x)/sinh(x);
  }
  else if (K==1){
    CotK = 1.0/tan(x);
    Phim2 = 1.0/beta*sin(beta*x)/sin(x);
  }
  else{
    CotK = 1.0/x;
    Phim2 = 1.0/beta*sin(beta*x)/x;
  }
  if (l==0){
    *Phi_l = Phim2;
    return _SUCCESS_;
  }
  Phim1 = Phim2*(CotK-beta/tan(beta*x))/sqrt(beta2-K);

 for (ll=2; ll<=l; ll++){
    Phi = ((2*ll-1)*CotK*Phim1-sqrt(beta2-K*(ll-1.0)*(ll-1.0))*Phim2)/
      sqrt(beta2-K*ll*ll);
    Phim2 = Phim1;
    Phim1 = Phi;
 }
 *Phi_l = Phim1; //Phim1 instead of Phi to take care of l==1 case.

 return _SUCCESS_;
}

int HypersphericalBackward(int K,int l,double beta,double x, double *Phi_l){
 double beta2 = beta*beta;
 double Phi0, Phi, Phipr1, Phim1,Phip1;
 double CotK;
 int isign, ll;
 double fj;

  if (K==-1){
    CotK = 1.0/tanh(x);
    Phi0 = 1.0/beta*sin(beta*x)/sinh(x);
  }
  else if (K==1){
    CotK = 1.0/tan(x);
    Phi0 = 1.0/beta*sin(beta*x)/sin(x);
  }
  else{
    CotK = 1.0/x;
    Phi0 = 1.0/beta*sin(beta*x)/x;
  }

  get_CF1(K, l, beta, x, &fj, &isign);

  *Phi_l = isign;
  Phipr1 = fj*(*Phi_l);
  Phi = (*Phi_l);
  Phip1 = 1.0/sqrt(beta2-K*(l+1.0)*(l+1.0))*(l*CotK*Phi-Phipr1);

  for (ll=l; ll>0; ll--){
    Phim1 = ( (2*ll+1)*CotK*Phi-sqrt(beta2-K*(ll+1.0)*(ll+1.0))*Phip1 )/
      sqrt(beta2-K*ll*ll);
    Phip1 = Phi;
    Phi = Phim1;
  }

  *Phi_l *= Phi0/Phi;

  return _SUCCESS_;
}



double airy_cheb_approx(double z){
  double Ai;
  if (z<=-7){
    Ai = coef1(z);
    return Ai;
  }
  if (z<=0){
    Ai = coef2(z);
    return Ai;
  }
  if (z<7){
    Ai = coef3(z);
  return Ai;
  }
  Ai = coef4(z);
  return Ai;
}

double coef1(double z){
  double A[5] = {1.1282427601,-0.6803534e-4,0.16687e-6,-0.128e-8,0.2e-10};
  double B[5] = {0.7822108673e-1,-0.6895649e-4,0.32857e-6,-0.37e-8,0.7e-10};
  double x,y,t,Ai,zeta,theta,sintheta,costheta,FA,FB;
 
  x = -z;
  zeta = 2.0/3.0*x*sqrt(x);
  theta = zeta+0.25*_PI_;
  sintheta = sin(theta);
  costheta = cos(theta);

  y = (7.0/x)*(7.0/x)*(7.0/x);
  FA = cheb(y,5,A);
  FB = cheb(y,5,B)/zeta;

  t = pow(x,-0.25);
  Ai = t*(sintheta*FA-costheta*FB);
  //Bi = t*(costheta*FA+sintheta*FB);
  return Ai;
}

double coef2(double z){
  double A[17] = {0.11535880704,0.6542816649e-1,0.26091774326,0.21959346500,
	      0.12476382168,-0.43476292594,0.28357718605,-0.9751797082e-1,
	      0.2182551823e-1,-0.350454097e-2,0.42778312e-3,
	      -0.4127564e-4,0.323880e-5,-0.21123e-6,0.1165e-7,
	      -0.55e-9,0.2e-10};
  double B[16] = {0.10888288487,-0.17511655051,0.13887871101,-0.11469998998,
	     0.22377807641,-0.18546243714,0.8063565186e-1,
	     -0.2208847864e-1,0.422444527e-2,-0.60131028e-3,
	     0.6653890e-4,-0.590842e-5,0.43131e-6,-0.2638e-7,
	     0.137e-8,-0.6e-10};
  //Ej = 3^(-j/3)/Gamma(j/3);
  double E1 = 0.355028053887817, E2 = 0.258819403792807;
  double x,FA,FB,Ai;
  x = -(z/7.0)*(z/7.0)*(z/7.0);
  FA = E1*cheb(x,17,&(A[0]));
  FB = E2*z*cheb(x,16,&(B[0]));
  Ai = FA-FB;
  //Bi = sqrt(3)*(FA+FB);
  return Ai;
}
double coef3(double z){
  double A[20] = {1.2974695794,-0.20230907821,
	      -0.45786169516,0.12953331987,0.6983827954e-1,
	      -0.3005148746e-1,-0.493036334e-2,0.390425474e-2,
	      -0.1546539e-4,-0.32193999e-3,0.3812734e-4,0.1714935e-4,
	      -0.416096e-5,-0.50623e-6,0.26346e-6,-0.281e-8,
	      -0.1122e-7,0.120e-8,0.31e-9,-0.7e-10};
  /*
  double B[25]={0.47839902387,-0.6881732880e-1,0.20938146768,
	    -0.3988095886e-1,0.4758441683e-1,-0.812296149e-2,
	    0.462845913e-2,0.70010098e-3,-0.75611274e-3,
	    0.68958657e-3,-0.33621865e-3,0.14501668e-3,-0.4766359e-4,
	    0.1251965e-4,-0.193012e-5,-0.19032e-6,0.29390e-6,
	    -0.13436e-6,0.4231e-7,-0.967e-8,0.135e-8,0.7e-10,
	    -0.12e-9,0.4e-10,-0.1e-10};
  */
  double x,EX,EY,Ai;
  x = z/7.0;
  EX = exp(1.75*z);
  EY = 1.0/EX;

  Ai = EY*cheb(x,20,&(A[0]));
  //Bi = EX*cheb(x,25,&(B[0]));
  return Ai;
}

double coef4(double z){
  double A[7]={0.56265126169,-0.76136219e-3,0.765252e-5,-0.14228e-6,
	    0.380e-8,-0.13e-9,0.1e-10};
  /*
  double B[7]={1.1316635302,0.166141673e-02,0.1968882e-04,0.47047e-06,
	    0.1769e-7,0.94e-9,0.6e-10};
  */
  double x,Y,t,zeta,EX,EY,Ai;

  Y = z*sqrt(z);
  zeta = 2.0/3.0*Y;
  EX = exp(zeta);
  EY = 1.0/EX;
  x = 7*sqrt(7)/Y;
  t = pow(z,-0.25);

  Ai = t*EY*cheb(x, 7, &(A[0]));
  //Bi = t*EX*cheb(x, 7, &(B[0]));
  return Ai;
}

double cheb(double x, int n, double *A){
  double b,d,u,y,c,F;
  int j;
  b = 0.0;
  d = A[n-1];
  u = 2*x-1.0;
  y = 2*u;
  
  for (j=(n-1);j>1;j--){
    c=b;
    b=d;
    d=y*b-c+A[j-1];
  }
  F = u*d-b+0.5*A[0];
  return F;
}


double get_value_at_small_phi(int K,int l,double beta,double Phi){
  double nu, lhs, alpha, xval;

  nu = l+0.5;
  lhs = 1.0/nu*log(2*Phi*nu);
  alpha = -2*lhs/5.0*(1.0+2.0*cosh(1.0/3.0*acosh(1.0+375/(16.0*lhs*lhs))));
  xval = nu/cosh(alpha)/beta;
  //Correct for geometry:
  if (K==1)
    xval *= asin(l/beta)/(l/beta);
  else if(K==-1)
    xval *= asinh(l/beta)/(l/beta);
  return xval;
}

int ClosedModY(int l, int beta, double *y){
  int phisign = 1;
  *y = fmod(*y,2.0*_PI_);
  if ((*y) > _PI_){
    *y = 2.0*_PI_-(*y);
    //phisign *= pow(-1,l)
    if (l%2==1) //l is uneven
      phisign = -phisign;
  }
  if ((*y)>0.5*_PI_){
    *y = _PI_-(*y);
    //phisign *= pow(-1,beta-l-1)
    if ((beta-l)%2==0) //beta-l-1 uneven
      phisign = -phisign;
  }
  return phisign;
}
