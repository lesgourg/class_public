/** @file hyperspherical.c Documented hyperspherical bessel function module.
 *
 * Thomas Tram, 11.01.2013
 *
 * This module computes hyperspherical Bessel functions.
 */

#include "hyperspherical.h"

int hyperspherical_HIS_create(int K, 
                              double beta, 
                              int nl, 
                              int *lvec, 
                              double xmin, 
                              double xmax, 
                              double sampling,
                              HyperInterpStruct **ppHIS, 
                              ErrorMsg error_message){
  /** Allocate storage for Hyperspherical Interpolation Structure (HIS).
      A pointer to a HIS pointer is input, so the structure itself is 
      also allocated here.
      Then, compute the values of Phi and dPhi and complete the interpolation
      structure.
  */
  HyperInterpStruct *pHIS;
  double deltax, beta2, lambda, x, xfwd, relerr;
  double *sqrtK, *PhiL;
  int j, k, l, nx, lmax;
  beta2 = beta*beta;
  lmax = lvec[nl-1];
  /*
  xmin = max(xmin,_HYPER_SAFETY_);
  if (K==1) 
    xmax = min(xmax,_PI_/2.0-_HYPER_SAFETY_); //We only need solution on [0;pi/2]
  */
  lambda = 2*_PI_/(beta+5.0); //Just to prevent too sparse sampling at beta<5.
  nx = (int) ((xmax-xmin)*sampling/lambda);
  nx = max(nx,2);
  deltax = (xmax-xmin)/(nx-1.0);
  //fprintf(stderr,"dx=%e\n",deltax);
  // Allocate vectors:
  class_alloc(pHIS,sizeof(HyperInterpStruct),error_message);
  class_alloc(pHIS->l,sizeof(int)*nl,error_message);
  class_alloc(pHIS->x,sizeof(double)*nx,error_message);
  class_alloc(pHIS->sinK,sizeof(double)*nx,error_message);
  class_alloc(pHIS->cotK,sizeof(double)*nx,error_message);
  class_alloc(pHIS->phi,sizeof(double)*nx*nl,error_message);
  class_alloc(pHIS->dphi,sizeof(double)*nx*nl,error_message);
  //Assign pHIS pointervalue to input ppHIS:
  *ppHIS = pHIS;
  //Set scalar values:
  pHIS->beta = beta;
  pHIS->delta_x = deltax;
  pHIS->l_size = nl;
  pHIS->x_size = nx;
  pHIS->K = K;
  //Order needed for trig interpolation: (We are using Taylor's remainder theorem)
  if (0.5*deltax*deltax < _TRIG_PRECISSION_)
    pHIS->trig_order = 1;
  else if ((pow(deltax,4)/24.0) < _TRIG_PRECISSION_)
    pHIS->trig_order = 3;
  else
    pHIS->trig_order = 5;

  //Copy lvector:
  for (j=0; j<nl; j++){
    pHIS->l[j] = lvec[j];
  }
  //Allocate sqrtK, and PhiL:
  class_alloc(sqrtK,(lmax+3)*sizeof(double),error_message);
  class_alloc(PhiL,(lmax+2)*sizeof(double),error_message);
  
  //Create xvector and set x, cotK, sinK, sqrtK and fwdidx:
  if (K==0){
    xfwd = sqrt(lmax*(lmax+1.0))/beta;
    for (j=0; j<nx; j++){
      x = xmin + j*deltax;
      pHIS->x[j] = x;
      pHIS->sinK[j] = x;
      pHIS->cotK[j] = 1.0/x;
    }
    for (l=0; l<=(lmax+2); l++){
      sqrtK[l] = beta;
    }
  }
  else if (K==1){
    xfwd = xmax+1.0;
    for (j=0; j<nx; j++){
      x = xmin + j*deltax;
      pHIS->x[j] = x;
      pHIS->sinK[j] = sin(x);
      pHIS->cotK[j] = 1.0/tan(x);
    }
    for (l=0; l<=(lmax+2); l++){
      sqrtK[l] = sqrt(beta2-l*l);
    }
    if (((int) (beta+0.2))==(lmax+1)){
      /** Take care of special case lmax = beta-1. The routine below will try to compute
          Phi_{lmax+1} which is not allowed. However, the purpose is to calculate the derivative
          Phi'_{lmax}, and the formula is correct if we set Phi_{lmax+1} = 0.
      */
      PhiL[lmax+1] = 0.0;
      lmax--;
    }
  }
  else if (K==-1){
    xfwd = asinh(sqrt(lmax*(lmax+1.0))/beta);
    for (j=0; j<nx; j++){
      x = xmin + j*deltax;
      pHIS->x[j] = x;
      pHIS->sinK[j] = sinh(x);
      pHIS->cotK[j] = 1.0/tanh(x);
    }
    for (l=0; l<=(lmax+2); l++){
      sqrtK[l] = sqrt(beta2+l*l);
    }
  }  
  //Calculate and assign Phi and dPhi values:
  for (j=0; j<nx; j++){
    x = pHIS->x[j];
    if (x<xfwd){
      //Use backwards method:
      hyperspherical_backwards_recurrence(K, 
                                          lmax+1, 
                                          beta, 
                                          x, 
                                          pHIS->sinK[j],
                                          pHIS->cotK[j],
                                          sqrtK,
                                          PhiL);
    }
    else{
      //Use forwards method:
      hyperspherical_forwards_recurrence(K, 
                                         lmax+1, 
                                         beta, 
                                         x, 
                                         pHIS->sinK[j],
                                         pHIS->cotK[j],
                                         sqrtK,
                                         PhiL);
    }
    //We have now populated PhiL at x, assign Phi and dPhi for all l in lvec:
    for (k=0; k<nl; k++){
      l = lvec[k];
      pHIS->phi[k*nx+j] = PhiL[l];
      pHIS->dphi[k*nx+j] = l*pHIS->cotK[j]*PhiL[l]-sqrtK[l+1]*PhiL[l+1];
      //      printf("x = %g, Phi_%d = %g %g\n",x,l,pHIS->phi[k*nx+j],pHIS->dphi[k*nx+j]);
    }
  }
  free(sqrtK);
  free(PhiL);
  return _SUCCESS_;
}


int hyperspherical_HIS_free(HyperInterpStruct *pHIS){
  /** Free the Hyperspherical Interpolation Structure. */  
  free(pHIS->l);
  free(pHIS->x);
  free(pHIS->sinK);
  free(pHIS->cotK);
  free(pHIS->phi);
  free(pHIS->dphi);
  free(pHIS);
  return _SUCCESS_;
}

int hyperspherical_Hermite_interpolation_vector(HyperInterpStruct *pHIS,
                                                int nxi,
                                                int lnum,
                                                double *xinterp,
                                                double *Phi,
                                                double *dPhi,
                                                double *d2Phi,
                                                double *sinKinterp,
                                                double *cosKinterp){
  /** Hermite interpolation of order 6 for Phi, dPhi, and d2Phi. When xinterp
      is sorted (increasing), computations can be reused. On the other hand, 
      for a randomly called value, the routine is not much slower than a
      routine optimised for this case. The more sorted the vector, the faster
      the execution time. For closed case, the interpolation structure only 
      covers [safety;pi/2-safety]. The calling routine should respect this.
      if sinK and cosK are not NULL, we will also interpolate them.
  */

  int do_function=_TRUE_, do_first_derivative=_TRUE_;
  int do_second_derivative=_TRUE_, do_first_or_second_derivative=_TRUE_;
  int do_trig_linear = _FALSE_, do_trig_three = _FALSE_, do_trig_five = _FALSE_;
  double ym=0, yp=0, dym=0, dyp=0, d2ym=0, d2yp=0, x, z, z2, z3, z4, z5;
  double cotKm=0,cotKp=0,sinKm=0,sinKp=0, sinKm2, sinKp2;
  double d3ym = 0, d3yp=0, d4ym=0, d4yp=0;
  double a1=0, a2=0, a3=0, a4=0, a5=0;
  double b1=0, b2=0, b3=0, b4=0, b5=0;
  double c1=0, c2=0, c3=0, c4=0, c5=0;
  double d1=0, d2=0, d3=0, d4=0, d5=0;
  double e1=0, e2=0, e3=0, e4=0, e5=0;
  double beta, beta2, *xvec, *sinK, *cotK;
  double xmin, xmax, deltax, deltax2, lxlp1;
  double left_border, right_border, next_border;
  int K, l, j, nx, current_border_idx;
  double *Phi_l, *dPhi_l;
  
  /** Set logical flags. The compiler should probably generate 2^3-1=7
      different functions, according to these flags. If not, maybe I should
      do it.
  */
 
  
  if (Phi==NULL)
    do_function = _FALSE_;
  else
    do_function = _TRUE_;
  if (dPhi == NULL)
    do_first_derivative = _FALSE_;
  else
    do_first_derivative = _TRUE_;
  if (d2Phi == NULL)
    do_second_derivative = _FALSE_;
  else
    do_second_derivative = _TRUE_;
  if ((do_first_derivative == _TRUE_)||(do_second_derivative == _TRUE_))
    do_first_or_second_derivative = _TRUE_;
  else
    do_first_or_second_derivative = _FALSE_;
  if ((sinKinterp!=NULL)&&(cosKinterp!=NULL)){
    //Select interpolation method:
    if (pHIS->trig_order==5)
      do_trig_five = _TRUE_;
    else if(pHIS->trig_order==3)
      do_trig_three = _TRUE_;
    else
      do_trig_linear = _TRUE_;
  }

  xvec = pHIS->x;
  sinK = pHIS->sinK;
  cotK = pHIS->cotK;
  beta = pHIS->beta;
  beta2 = beta*beta;
  deltax = pHIS->delta_x;
  deltax2 = deltax*deltax;
  K = pHIS->K;
  nx = pHIS->x_size;
  Phi_l = pHIS->phi+lnum*nx;
  dPhi_l = pHIS->dphi+lnum*nx;
  l = pHIS->l[lnum];
  lxlp1 = l*(l+1.0);
  xmin = xvec[0];
  xmax = xvec[nx-1];

  left_border = xmax;
  right_border = xmin;
  next_border = xmin;

  for (j=0; j<nxi; j++){
    x = xinterp[j];
    //Loop over output values
    if ((x<xmin)||(x>xmax)){
      //Outside interpolation region, set to zero.
      if (do_function==_TRUE_) 
        Phi[j] = 0.0;
      if (do_first_derivative==_TRUE_)
        dPhi[j] = 0.0;
      if (do_second_derivative==_TRUE_)
        d2Phi[j] = 0.0;
      continue;
    }
    if ((x>right_border)||(x<left_border)){
      if ((x>next_border)||(x<left_border)){
        current_border_idx = ((int) ((x-xmin)/deltax))+1; 
        current_border_idx = max(1,current_border_idx);
        current_border_idx = min(nx-1,current_border_idx);
        //printf("Current border index at jump: %d\n",current_border_idx);
        //max operation takes care of case x = xmin,
        //min operation takes care of case x = xmax.
        //Calculate left derivatives:
        cotKm = cotK[current_border_idx-1];
        sinKm = sinK[current_border_idx-1];
        sinKm2 = sinKm*sinKm;
        ym = Phi_l[current_border_idx-1];
        dym = dPhi_l[current_border_idx-1];
        d2ym = -2*dym*cotKm+ym*(lxlp1/sinKm2-beta2+K);
        //printf("%g %g %g %g %g\n",cotKm,sinKm,ym,dym,d2ym);
        if (do_first_or_second_derivative==_TRUE_){
          d3ym = -2*cotKm*d2ym-2*ym*lxlp1*cotKm/sinKm2+
            dym*(K-beta2+(2+lxlp1)/sinKm2);
        }
        if (do_second_derivative==_TRUE_){
          d4ym = -2*cotKm*d3ym + d2ym*(K-beta2+(4+lxlp1)/sinKm2)+
            dym*(-4*(1+lxlp1)*cotKm/sinKm2)+
            ym*(2*lxlp1/sinKm2*(2*cotKm*cotKm+1/sinKm2));
        }
      }
      else{
        //x>current_border but not next border: I have moved to next block.
        current_border_idx++;
        //printf("Current border index at else: %d\n",current_border_idx);
        //Copy former right derivatives to left derivatives.
        ym = yp;
        dym = dyp;
        d2ym = d2yp;
        d3ym = d3yp;
        d4ym = d4yp;
        sinKm = sinKp;
        cotKm = cotKp;
      }
      left_border = xvec[max(0,current_border_idx-1)];
      right_border = xvec[current_border_idx];
      next_border = xvec[min(nx-1,current_border_idx+1)];
      //Evaluate right derivatives and calculate coefficients:
      cotKp = cotK[current_border_idx];
      sinKp = sinK[current_border_idx];
      sinKp2 = sinKp*sinKp;
      yp = Phi_l[current_border_idx];
      dyp = dPhi_l[current_border_idx];
      d2yp = -2*dyp*cotKp+yp*(lxlp1/sinKp2-beta2+K);
      if (do_first_or_second_derivative == _TRUE_){
        d3yp = -2*cotKp*d2yp-2*yp*lxlp1*cotKp/sinKp2+
          dyp*(K-beta2+(2+lxlp1)/sinKp2);
      }
      if (do_second_derivative == _TRUE_){
        d4yp = -2*cotKp*d3yp + d2yp*(K-beta2+(4+lxlp1)/sinKp2)+
          dyp*(-4*(1+lxlp1)*cotKp/sinKp2)+
          yp*(2*lxlp1/sinKp2*(2*cotKp*cotKp+1/sinKp2));
      }
      if (do_function == _TRUE_){
        a1 = dym*deltax;
        a2 = 0.5*d2ym*deltax2;
        a3 = (-1.5*d2ym+0.5*d2yp)*deltax2-(6*dym+4*dyp)*deltax-10*(ym-yp);
        a4 = (1.5*d2ym-d2yp)*deltax2+(8*dym+7*dyp)*deltax+15*(ym-yp);
        a5 = (-0.5*d2ym+0.5*d2yp)*deltax2-3*(dym+dyp)*deltax-6*(ym-yp);
      }
      if (do_first_derivative==_TRUE_){
        b1 = d2ym*deltax;
        b2 = 0.5*d3ym*deltax2;
        b3 = (-1.5*d3ym+0.5*d3yp)*deltax2-(6*d2ym+4*d2yp)*deltax-10*(dym-dyp);
        b4 = (1.5*d3ym-d3yp)*deltax2+(8*d2ym+7*d2yp)*deltax+15*(dym-dyp);
        b5 = (-0.5*d3ym+0.5*d3yp)*deltax2-3*(d2ym+d2yp)*deltax-6*(dym-dyp);
      }
      if (do_second_derivative==_TRUE_){
        c1 = d3ym*deltax;
        c2 = 0.5*d4ym*deltax2;
        c3 = (-1.5*d4ym+0.5*d4yp)*deltax2-(6*d3ym+4*d3yp)*deltax-10*(d2ym-d2yp);
        c4 = (1.5*d4ym-d4yp)*deltax2+(8*d3ym+7*d3yp)*deltax+15*(d2ym-d2yp);
        c5 = (-0.5*d4ym+0.5*d4yp)*deltax2-3*(d3ym+d3yp)*deltax-6*(d2ym-d2yp);
      }
      if (do_trig_five==_TRUE_){
        //sinK
        d1 = sinKm*cotKm*deltax;
        d2 = -K*0.5*sinKm*deltax2;
        d3 = -K*(-1.5*sinKm+0.5*sinKp)*deltax2-(6*sinKm*cotKm+4*sinKp*cotKp)*deltax
          -10*(sinKm-sinKp);
        d4 = -K*(1.5*sinKm-sinKp)*deltax2+(8*sinKm*cotKm+7*sinKp*cotKp)*deltax+15*(sinKm-sinKp);
        d5 = -K*(-0.5*sinKm+0.5*sinKp)*deltax2-3*(sinKm*cotKm+sinKp*cotKp)*deltax-6*(sinKm-sinKp);
        //cosK:
        e1 = -K*sinKm*deltax;
        e2 = -K*0.5*sinKm*cotKm*deltax2;
        e3 = -K*(-1.5*sinKm*cotKm+0.5*sinKp*cotKp)*deltax2+K*(6*sinKm+4*sinKp)*deltax-
          10*(cotKm*sinKm-cotKp*sinKp);
        e4 = -K*(1.5*sinKm*cotKm-sinKp*cotKp)*deltax2-K*(8*sinKm+7*sinKp)*deltax+
          15*(cotKm*sinKm-cotKp*sinKp);
        e5 = -K*(-0.5*sinKm*cotKm+0.5*sinKp*cotKp)*deltax2+
          3*K*(sinKm+sinKp)*deltax-6*(cotKm*sinKm-cotKp*sinKp);
      }
      else if (do_trig_three==_TRUE_){
        //sin
        d1 = sinKm*cotKm*deltax;
        d2 = -2*sinKm*cotKm*deltax-sinKp*cotKp*deltax-3*sinKm+3*sinKp;
        d3 = sinKm*cotKm*deltax+sinKp*cotKp*deltax+2*sinKm-2*sinKp;
        //cos
        e1 = -K*sinKm*deltax;
        e2 = 2*K*sinKm*deltax+K*sinKp*deltax-3*sinKm*cotKm+3*sinKp*cotKp;
        e3 = -K*sinKm*deltax-K*sinKp*deltax+2*sinKm*cotKm-2*sinKp*cotKp;
      }
    }
    //Evaluate polynomial:
    z = (x-left_border)/deltax;
    z2 = z*z;
    z3 = z2*z;
    z4 = z2*z2;
    z5 = z2*z3;
    if (do_function == _TRUE_)
      Phi[j] = ym+a1*z+a2*z2+a3*z3+a4*z4+a5*z5;
    if (do_first_derivative == _TRUE_)
      dPhi[j] = dym+b1*z+b2*z2+b3*z3+b4*z4+b5*z5;
    if (do_second_derivative == _TRUE_)
      d2Phi[j] = d2ym+c1*z+c2*z2+c3*z3+c4*z4+c5*z5;
    if (do_trig_five==_TRUE_){
      sinKinterp[j] = sinKm+d1*z+d2*z2+d3*z3+d4*z4+d5*z5;
      cosKinterp[j] = cotKm*sinKm+e1*z+e2*z2+e3*z3+e4*z4+e5*z5;
    } else if (do_trig_three==_TRUE_){
      sinKinterp[j] = sinKm+d1*z+d2*z2+d3*z3;
      cosKinterp[j] = cotKm*sinKm+e1*z+e2*z2+e3*z3;
    } else if (do_trig_linear==_TRUE_){
      //We do linear interpolation directly, sinKm, sinKp are available
      sinKinterp[j] = sinKm+(sinKp-sinKm)*z;
      cosKinterp[j] = cotKm*sinKm*(1.0-z)+cotKp*sinKp*z;
    }
    //printf("x = %g, [%g, %g, %g]\n",x,Phi[j],dPhi[j],d2Phi[j]);
  }
  return _SUCCESS_;
}

int hyperspherical_forwards_recurrence(int K, 
                                       int lmax, 
                                       double beta, 
                                       double x, 
                                       double sinK,
                                       double cotK,
                                       double *sqrtK,
                                       double *PhiL){
  int l;

  //  printf("x = %g, K=%d, beta = %g, lmax = %d, sinK = %g, cotK = %g.\n",
  //     x,K,beta,lmax,sinK,cotK);
  PhiL[0] = 1.0/beta*sin(beta*x)/sinK;
  PhiL[1] = PhiL[0]*(cotK-beta/tan(beta*x))/sqrtK[1];
  for (l=2; l<=lmax; l++){
    PhiL[l] = ((2*l-1)*cotK*PhiL[l-1]-PhiL[l-2]*sqrtK[l-1])/sqrtK[l];
  }
  //printf("Phi_0 = %g, Phi_1=%g, Phi_lmax = %g\n",PhiL[0],PhiL[1], PhiL[lmax]);
  return _SUCCESS_;
}

int hyperspherical_backwards_recurrence(int K, 
                                        int lmax, 
                                        double beta, 
                                        double x, 
                                        double sinK,
                                        double cotK,
                                        double *sqrtK,
                                        double *PhiL){
  double phi0, phi1, phipr1, phi, phi_plus_1_times_sqrtK, phi_minus_1, scaling;
  int l, k, isign;
  phi0 = sin(beta*x)/(beta*sinK);

  //printf("in backwards. x = %g\n",x);
  if (K==1){
    CF1_from_Gegenbauer(lmax,(int) (beta+0.2),sinK,cotK, &phipr1);
    phi1 = 1.0;
  }
  else{
    get_CF1(K,lmax,beta,cotK, &phipr1, &isign);
    phi1 = isign;
    phipr1 *=phi1;  
    //printf("isign = %d, phi1 = %g, phipr1 = %g\n",isign,phi1,phipr1);
  }

  PhiL[lmax] = phi1;
  phi = phi1;
  //  phi_plus_1 = 1/sqrtK[lmax+1]*(lmax*cotK*phi1-phipr1);
  phi_plus_1_times_sqrtK = lmax*cotK*phi1-phipr1;
  for (l=lmax; l>=1; l--){
    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
    phi_plus_1_times_sqrtK = phi*sqrtK[l];
    phi = phi_minus_1;
    PhiL[l-1] = phi;
    if (fabs(phi)>_HYPER_OVERFLOW_){
      phi = phi/_HYPER_OVERFLOW_;
      phi_plus_1_times_sqrtK = phi_plus_1_times_sqrtK/_HYPER_OVERFLOW_;
      //Rescale whole Phi vector until this point:
      for (k=l-1; k<=lmax; k++)
        PhiL[k] /= _HYPER_OVERFLOW_;
    }
  }
  scaling = phi0/phi;
  for (k=0; k<=lmax; k++)
    PhiL[k] *= scaling;
  return _SUCCESS_;
}

int get_CF1(int K,int l,double beta, double cotK, double *CF, int *isign){
  int maxiter = 1000000;
  double tiny = 1e-100;
  double reltol = DBL_EPSILON;
  double aj,bj,fj,Cj,Dj,Delj;
  double beta2 = beta*beta;
  double sqrttmp;
  int j;

  bj = l*cotK; //This is b_0
  fj = bj;
  Cj = bj;
  Dj = 0.0;
  *isign = 1;
  for(j=1; j<=maxiter; j++){
    sqrttmp = sqrt(beta2-K*(l+j+1)*(l+j+1));
    aj = -sqrt(beta2-K*(l+j)*(l+j))/sqrttmp;
    if (j==1)
      aj = sqrt(beta2-K*(l+1)*(l+1))*aj;
    bj = (2*(l+j)+1)/sqrttmp*cotK;
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
      //printf("iter:%d, %g, %g\n",j,sqrttmp,fj);
      return _SUCCESS_;
    }
  }
  return _FAILURE_;
}

int CF1_from_Gegenbauer(int l,
                        int beta, 
                        double sinK, 
                        double cotK, 
                        double *CF){
  int n, alpha, k;
  double x, G, dG, Gkm1, Gkm2;
  n = beta-l-1;
  alpha = l+1;
  x = sinK*cotK; //cos(x)
  switch (n){
  case 0:
    G = 1;
    dG = 0;
    break;
  case 1:
    G = 2*alpha*x;
    dG = 2*alpha;
    break;
  case 2:
    G = -alpha + 2*alpha*(1+alpha)*x*x;
    dG = 4*x*alpha*(1+alpha);
    break;
  case 3:
    G = -2*alpha*(1+alpha)*x+4.0/3.0*alpha*(1+alpha)*(2+alpha)*x*x*x;
    dG = 2*alpha*(1+alpha)*(2*(2+alpha)*x*x-1);
    break;
  default:
    Gkm2 = -alpha + 2*alpha*(1+alpha)*x*x;
    Gkm1 = -2*alpha*(1+alpha)*x+4.0/3.0*alpha*(1+alpha)*(2+alpha)*x*x*x;
    for (k=4; k<=n; k++){
      G = (2*(k+alpha-1)*x*Gkm1 - (k+2*alpha-2)*Gkm2) / k;
      if (fabs(G)>_HYPER_OVERFLOW_){
                Gkm2 = Gkm1/_HYPER_OVERFLOW_;
                G = G/_HYPER_OVERFLOW_;
                Gkm1 = G;
      }
      else{
        Gkm2 = Gkm1;
        Gkm1 = G;
      }
    }
    //Derivative. Gkm2 is in fact Gkm1..
    dG = (-n*x*G+(n+2*alpha-1)*Gkm2)/(1.0-x*x);
  }
  //%Phi = G;
  //%dPhi = l*coty.*G-siny.*dG;
  *CF = l*cotK-sinK*dG/G;
}
 
int hyperspherical_WKB(int K,int l,double beta,double y, double *Phi){
  double e, w, w2, alpha, alpha2, CscK, ytp, t;
  double S, Q, C, argu, Ai;
  int airy_sign = 1, phisign = 1, intbeta;
  double ldbl = l;

  if (K==1){
    //Limit range to [0; pi/2]:
    intbeta = (int)(beta+0.4); //Round to nearest integer (just to be sure)
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
  double B[25]={0.47839902387,-0.6881732880e-1,0.20938146768,
            -0.3988095886e-1,0.4758441683e-1,-0.812296149e-2,
            0.462845913e-2,0.70010098e-3,-0.75611274e-3,
            0.68958657e-3,-0.33621865e-3,0.14501668e-3,-0.4766359e-4,
            0.1251965e-4,-0.193012e-5,-0.19032e-6,0.29390e-6,
            -0.13436e-6,0.4231e-7,-0.967e-8,0.135e-8,0.7e-10,
            -0.12e-9,0.4e-10,-0.1e-10};
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
  double B[7]={1.1316635302,0.166141673e-02,0.1968882e-04,0.47047e-06,
            0.1769e-7,0.94e-9,0.6e-10};
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


int hyperspherical_get_xmin(HyperInterpStruct *pHIS,
                            double xtol,
                            double phiminabs,
                            double *xmin){
  int left_index, right_index, index_l, j;
  int nl = pHIS->l_size;
  int nx = pHIS->x_size;
  int REFINE=10;
  double x[REFINE];
  double Phi[REFINE];
  double *phivec = pHIS->phi;
  double *xvec = pHIS->x;
  double xleft, xright;
  
  for (index_l=0; index_l<nl; index_l++){
    for (right_index = 0; right_index<nx; right_index++){
      if (fabs(phivec[index_l*nx+right_index])>phiminabs)
        break;
    }
    if (right_index==0){
      xmin[index_l] = xvec[0];
      //printf("special case: xmin = %.16e for index_l=%d\n",xmin[index_l],index_l);
      continue;
    }
    if (right_index==nx){
      xmin[index_l] = xvec[nx-1];
      //printf("special case: xmin = %.16e for index_l=%d\n",xmin[index_l],index_l);
      continue;
    }
    left_index = right_index-1;
    xleft = xvec[left_index];
    xright = xvec[right_index];
    xmin[index_l] = xright;
    while ((xright-xleft)>xtol){
      //Create interpolation vector
      //printf("Refining\n");
      for (j=0; j<REFINE; j++)
        x[j] = xleft+j*(xright-xleft)/(REFINE-1.0);
      hyperspherical_Hermite_interpolation_vector(pHIS,REFINE,
                                                  index_l, x, Phi, NULL,NULL,NULL,NULL);
      for (right_index = 1; right_index<REFINE; right_index++){
        if (fabs(Phi[right_index])>phiminabs)
          break;
      }
      left_index = right_index-1;
      xleft = x[left_index];
      xright = x[right_index];
      xmin[index_l] = xright;
    }
    //printf("xmin = %.16e\n",xmin[index_l]);
  }
  return _SUCCESS_;
}
