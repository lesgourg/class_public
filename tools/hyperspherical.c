/** @file hyperspherical.c Documented hyperspherical bessel function module.
 *
 * Thomas Tram, 11.01.2013
 *
 * This module computes hyperspherical Bessel functions.
 */

#include "hyperspherical.h"
#include "parallel.h"

int hyperspherical_HIS_create(int K,
                              double beta,
                              int nl,
                              int *lvec,
                              double xmin,
                              double xmax,
                              double sampling,
                              int l_WKB,
                              double phiminabs,
                              HyperInterpStruct *pHIS,
                              ErrorMsg error_message){
  /** Allocate storage for Hyperspherical Interpolation Structure (HIS).
      Then, compute the values of Phi and dPhi and complete the interpolation
      structure.
      Change to accomodate shared memory approach: Allocate all memory in a
      single call, and return the pointer as ppHIS. All pointers inside are
      then relative to ppHIS.
  */
  double deltax, beta2, lambda, x, xfwd;
  double *sqrtK, *one_over_sqrtK;
  int j, k, l, nx, lmax, l_recurrence_max;

  beta2 = beta*beta;
  lmax = lvec[nl-1];
  lambda = 2*_PI_/beta;
  nx = (int) ((xmax-xmin)*sampling/lambda);
  nx = MAX(nx,2);
  deltax = (xmax-xmin)/(nx-1.0);
  //fprintf(stderr,"dx=%e\n",deltax);
  //fprintf(stderr,"%e %e\n",beta,sampling);
  //Set scalar values:
  pHIS->beta = beta;
  pHIS->delta_x = deltax;
  pHIS->l_size = nl;
  pHIS->x_size = nx;
  pHIS->K = K;
  //Set pointervalues in pHIS:

  class_alloc(pHIS->l, sizeof(int)*nl,error_message);
  class_alloc(pHIS->chi_at_phimin,sizeof(double)*nl,error_message);
  class_alloc(pHIS->x,sizeof(double)*nx,error_message);
  class_alloc(pHIS->sinK,sizeof(double)*nx,error_message);
  class_alloc(pHIS->cotK,sizeof(double)*nx,error_message);
  class_alloc(pHIS->phi,sizeof(double)*nx*nl,error_message);
  class_alloc(pHIS->dphi,sizeof(double)*nx*nl,error_message);

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
  class_alloc(one_over_sqrtK,(lmax+3)*sizeof(double),error_message);
  //class_alloc(PhiL,(lmax+2)*sizeof(double),error_message);

  //Find l_WKB_min, the highest l in lvec where l<l_WKB:
  l_recurrence_max = -10;
  int index_recurrence_max=-10;
  for (k=nl-1; k>=0; k--){
    l = lvec[k];
    if (l<l_WKB){
        l_recurrence_max = l;
        index_recurrence_max = k;
        break;
    }
  }

  //Create xvector and set x, cotK, sinK, sqrtK and fwdidx:
  switch (K){
  case 0:
    xfwd = sqrt(l_recurrence_max*(l_recurrence_max+1.0))/beta;
    for (j=0; j<nx; j++){
      x = xmin + j*deltax;
      pHIS->x[j] = x;
      pHIS->sinK[j] = x;
      pHIS->cotK[j] = 1.0/x;
    }
    for (l=0; l<=(lmax+2); l++){
      sqrtK[l] = beta;
      one_over_sqrtK[l] = 1.0/sqrtK[l];
    }
    break;
  case 1:
    xfwd = asin(sqrt(l_recurrence_max*(l_recurrence_max+1.0))/beta);
    for (j=0; j<nx; j++){
      x = xmin + j*deltax;
      pHIS->x[j] = x;
      pHIS->sinK[j] = sin(x);
      pHIS->cotK[j] = 1.0/tan(x);
    }
    for (l=0; l<=(lmax+2); l++){
      sqrtK[l] = sqrt(beta2-l*l);
      one_over_sqrtK[l] = 1.0/sqrtK[l];
    }
    break;
  case -1:
    xfwd = asinh(sqrt(l_recurrence_max*(l_recurrence_max+1.0))/beta);
    for (j=0; j<nx; j++){
      x = xmin + j*deltax;
      pHIS->x[j] = x;
      pHIS->sinK[j] = sinh(x);
      pHIS->cotK[j] = 1.0/tanh(x);
    }
    for (l=0; l<=(lmax+2); l++){
      sqrtK[l] = sqrt(beta2+l*l);
      one_over_sqrtK[l] = 1.0/sqrtK[l];
    }
    break;
  default:
    return _FAILURE_;
  }

  int xfwdidx = (xfwd-xmin)/deltax;
  //Calculate and assign Phi and dPhi values:

  int TNUM_i;
  int TNUM_TOTAL;
  class_setup_parallel();

  /* Increase beyond threads of system to increase parallelization
   *  This will require more memory space to allocate the PhiL,
   *  but memory space is cheaper than computation time, usually.
   * Default increase factor: 6 */
  TNUM_TOTAL = task_system.get_num_threads() * 6;
  for(TNUM_i=0;TNUM_i<TNUM_TOTAL;++TNUM_i){
    class_run_parallel(=,
      int declare_list_of_variables_inside_parallel_region(j,k,l,index_x,current_chunk);
      double * PhiL;
      int lmax_local;
      lmax_local = lmax;
      class_alloc(PhiL,(lmax_local+2)*sizeof(double)*_HYPER_CHUNK_,error_message);

      if ((K == 1) && ((int)(beta+0.2) == (lmax_local+1))) {
        /** Take care of special case lmax_local = beta-1.
            The routine below will try to compute
            Phi_{lmax_local+1} which is not allowed. However,
            the purpose is to calculate the derivative
            Phi'_{lmax_local}, and the formula is correct if we set Phi_{lmax_local+1} = 0.
        */
        for(j=0; j<(lmax_local+2)*_HYPER_CHUNK_; j++){
          PhiL[j] = 0.0;
        }

        lmax_local--;
      }

      for (j=0; j<MIN(nx,xfwdidx); j++){
        if((j%TNUM_TOTAL)!=TNUM_i){continue;} //Only do the work assigned to TNUM_i
        //Use backwards method:
        hyperspherical_backwards_recurrence(K,
                                            MIN(l_recurrence_max,lmax_local)+1,
                                            beta,
                                            pHIS->x[j],
                                            pHIS->sinK[j],
                                            pHIS->cotK[j],
                                            sqrtK,
                                            one_over_sqrtK,
                                            PhiL);
        //We have now populated PhiL at x, assign Phi and dPhi for all l in lvec:
        for (k=0; k<=index_recurrence_max; k++){
          l = lvec[k];
          pHIS->phi[k*nx+j] = PhiL[l];
          pHIS->dphi[k*nx+j] = l*pHIS->cotK[j]*PhiL[l]-sqrtK[l+1]*PhiL[l+1];
        }
      }
      /**
      for (j=0; j<MIN(nx,xfwdidx); j+= _HYPER_CHUNK_){
        current_chunk = MIN(_HYPER_CHUNK_,MIN(nx,xfwdidx)-j);
        //Use backwards method:
        hyperspherical_backwards_recurrence_chunk(K,
                                                  MIN(l_recurrence_max,lmax_local)+1,
                                                  beta,
                                                  pHIS->x+j,
                                                  pHIS->sinK+j,
                                                  pHIS->cotK+j,
                                                  current_chunk,
                                                  sqrtK,
                                                  one_over_sqrtK,
                                                  PhiL);
        //We have now populated PhiL at x, assign Phi and dPhi for all l in lvec:
        for (k=0; k<=index_recurrence_max; k++){
          l = lvec[k];
          for (index_x=0; index_x<current_chunk; index_x++){
            pHIS->phi[k*nx+j+index_x] = PhiL[l*current_chunk+index_x];
            pHIS->dphi[k*nx+j+index_x] = l*pHIS->cotK[j+index_x]*
              PhiL[l*current_chunk+index_x]-
              sqrtK[l+1]*PhiL[(l+1)*current_chunk+index_x];
          }
        }
      }

      */

      for (j=xfwdidx; j<nx; j+=_HYPER_CHUNK_){
        if((((j-xfwdidx)/_HYPER_CHUNK_)%TNUM_TOTAL)!=TNUM_i){continue;} //Only do the work assigned to TNUM_i
        //Use forwards method:
        current_chunk = MIN(_HYPER_CHUNK_,nx-j);
        hyperspherical_forwards_recurrence_chunk(K,
                                                 MIN(l_recurrence_max,lmax_local)+1,
                                                 beta,
                                                 pHIS->x+j,
                                                 pHIS->sinK+j,
                                                 pHIS->cotK+j,
                                                 current_chunk,
                                                 sqrtK,
                                                 one_over_sqrtK,
                                                 PhiL);

        //We have now populated PhiL at x, assign Phi and dPhi for all l in lvec:
        for (k=0; k<=index_recurrence_max; k++){
          l = lvec[k];
          for (index_x=0; index_x<current_chunk; index_x++){
            pHIS->phi[k*nx+j+index_x] = PhiL[l*current_chunk+index_x];
            pHIS->dphi[k*nx+j+index_x] = l*pHIS->cotK[j+index_x]*
              PhiL[l*current_chunk+index_x]-
              sqrtK[l+1]*PhiL[(l+1)*current_chunk+index_x];
          }
        }
      }
      free(PhiL);
      return _SUCCESS_;
    );
  }
  class_finish_parallel();

  free(sqrtK);
  free(one_over_sqrtK);

  for (k=0; k<nl; k++){
    hyperspherical_get_xmin_from_approx(K,lvec[k],beta,0.,phiminabs,pHIS->chi_at_phimin+k,NULL);
  }

  //hyperspherical_get_xmin(pHIS,1.e-4,phiminabs,pHIS->chi_at_phimin);

  return _SUCCESS_;
}

size_t hyperspherical_HIS_size(int nl, int nx){
  return(sizeof(int)*nl+sizeof(double)*nl+3*sizeof(double)*nx+2*sizeof(double)*nx*nl);
}

int hyperspherical_update_pointers(HyperInterpStruct *pHIS_local,
                                   void * HIS_storage_shared){
  /** Assign pointers in pHIS: (Remember that pointer incrementation moves
      the number of bytes taken up by 1 variable of the type that the
      pointer points to. */
  int nx=pHIS_local->x_size;
  int nl=pHIS_local->l_size;
  pHIS_local->l = (int *) (HIS_storage_shared);
  pHIS_local->chi_at_phimin = (double *) (pHIS_local->l+nl);
  pHIS_local->x = pHIS_local->chi_at_phimin+nl;
  pHIS_local->sinK = pHIS_local->x + nx;
  pHIS_local->cotK = pHIS_local->sinK + nx;
  pHIS_local->phi = pHIS_local->cotK +nx;
  pHIS_local->dphi = pHIS_local->phi+nx*nl;

  return _SUCCESS_;
}

int hyperspherical_HIS_free(HyperInterpStruct *pHIS,
                            ErrorMsg error_message){
  /** Free the Hyperspherical Interpolation Structure. */
  free(pHIS->l);
  free(pHIS->chi_at_phimin);
  free(pHIS->x);
  free(pHIS->sinK);
  free(pHIS->cotK);
  free(pHIS->phi);
  free(pHIS->dphi);

  return _SUCCESS_;
}

int hyperspherical_Hermite_interpolation_vector(HyperInterpStruct *pHIS,
                                                int nxi,
                                                int lnum,
                                                double *xinterp,
                                                double *Phi,
                                                double *dPhi,
                                                double *d2Phi) {

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
  double ym=0, yp=0, dym=0, dyp=0, d2ym=0, d2yp=0, x, z, z2, z3, z4, z5;
  double cotKm=0,cotKp=0,sinKm=0,sinKp=0, sinKm2, sinKp2;
  double d3ym = 0, d3yp=0, d4ym=0, d4yp=0;
  double a1=0, a2=0, a3=0, a4=0, a5=0;
  double b1=0, b2=0, b3=0, b4=0, b5=0;
  double c1=0, c2=0, c3=0, c4=0, c5=0;
  double beta, beta2, *xvec, *sinK, *cotK;
  double xmin, xmax, deltax, deltax2, lxlp1;
  double left_border, right_border, next_border;
  int K, l, j, nx, current_border_idx=0;
  double *Phi_l, *dPhi_l;
  int phisign = 1, dphisign = 1;

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
    //take advantage of periodicity of functions in closed case
    if (pHIS->K==1)
      ClosedModY(pHIS->l[lnum], (int)(pHIS->beta+0.2), &x, &phisign, &dphisign);
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
        current_border_idx = MAX(1,current_border_idx);
        current_border_idx = MIN(nx-1,current_border_idx);
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
      left_border = xvec[MAX(0,current_border_idx-1)];
      right_border = xvec[current_border_idx];
      next_border = xvec[MIN(nx-1,current_border_idx+1)];
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
    }
    //Evaluate polynomial:
    z = (x-left_border)/deltax;
    z2 = z*z;
    z3 = z2*z;
    z4 = z2*z2;
    z5 = z2*z3;
    if (do_function == _TRUE_)
      Phi[j] = (ym+a1*z+a2*z2+a3*z3+a4*z4+a5*z5)*phisign;
    if (do_first_derivative == _TRUE_)
      dPhi[j] = (dym+b1*z+b2*z2+b3*z3+b4*z4+b5*z5)*dphisign;
    if (do_second_derivative == _TRUE_)
      d2Phi[j] = (d2ym+c1*z+c2*z2+c3*z3+c4*z4+c5*z5)*phisign;
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
                                       double * __restrict__ sqrtK,
                                       double * __restrict__ one_over_sqrtK,
                                       double * __restrict__ PhiL){
  int l;
  PhiL[0] = 1.0/beta*sin(beta*x)/sinK;
  PhiL[1] = PhiL[0]*(cotK-beta/tan(beta*x))*one_over_sqrtK[1];
  for (l=2; l<=lmax; l++){
    PhiL[l] = ((2*l-1)*cotK*PhiL[l-1]-PhiL[l-2]*sqrtK[l-1])*one_over_sqrtK[l];
  }
  return _SUCCESS_;
}

int hyperspherical_forwards_recurrence_chunk(int K,
                                             int lmax,
                                             double beta,
                                             double * __restrict__ x,
                                             double * __restrict__ sinK,
                                             double * __restrict__ cotK,
                                             int chunk,
                                             double * __restrict__ sqrtK,
                                             double * __restrict__ one_over_sqrtK,
                                             double * __restrict__ PhiL){
  int l;
  int index_x;
  for (index_x=0; index_x<chunk; index_x++){
    PhiL[index_x] = 1.0/beta*sin(beta*x[index_x])/sinK[index_x];
    PhiL[chunk+index_x] = PhiL[index_x]*
      (cotK[index_x]-beta/tan(beta*x[index_x]))*one_over_sqrtK[1];
  }
  for (l=2; l<=lmax; l++){
    for (index_x=0; index_x<chunk; index_x++)
      PhiL[l*chunk+index_x] =
        ((2*l-1)*cotK[index_x]*PhiL[(l-1)*chunk+index_x]-
         PhiL[(l-2)*chunk+index_x]*sqrtK[l-1])*one_over_sqrtK[l];
  }
  return _SUCCESS_;
}


int hyperspherical_backwards_recurrence(int K,
                                        int lmax,
                                        double beta,
                                        double x,
                                        double sinK,
                                        double cotK,
                                        double * __restrict__ sqrtK,
                                        double * __restrict__ one_over_sqrtK,
                                        double * __restrict__ PhiL){
  double phi0, phi1, phipr1, phi, phi_plus_1_times_sqrtK, phi_minus_1, scaling;
  int l, k, isign;
  int funcreturn = _FAILURE_;
  phi0 = sin(beta*x)/(beta*sinK);

  //printf("in backwards. x = %g\n",x);
  if (K==1){
    if (beta > 1.5*lmax) {
      funcreturn = get_CF1(K,lmax,beta,cotK, &phipr1, &isign);
    }
    if (funcreturn == _FAILURE_) {
      CF1_from_Gegenbauer(lmax,(int) (beta+0.2),sinK,cotK, &phipr1);
    }
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


  int l_ini, l_align;
  l_align = lmax-lmax%_HYPER_BLOCK_;

  // Bring l down to _HYPER_BLOCK_ aligned region:
  for (l=lmax; l>l_align; l--){
    //    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )*one_over_sqrtK[l];
    phi_plus_1_times_sqrtK = phi*sqrtK[l];
    phi = phi_minus_1;
    PhiL[l-1] = phi;
  }
  for (l_ini=l_align; l_ini>0; l_ini -= _HYPER_BLOCK_){
    for (l=l_ini; l>(l_ini-_HYPER_BLOCK_); l--){
      //    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
      phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )*one_over_sqrtK[l];
      phi_plus_1_times_sqrtK = phi*sqrtK[l];
      phi = phi_minus_1;
      PhiL[l-1] = phi;
    }
    if (fabs(phi)>_HYPER_OVERFLOW_){
      phi *= _ONE_OVER_HYPER_OVERFLOW_;
      phi_plus_1_times_sqrtK *= _ONE_OVER_HYPER_OVERFLOW_;
      //Rescale whole Phi vector until this point:
      for (k=l; k<=lmax; k++)
        PhiL[k] *=_ONE_OVER_HYPER_OVERFLOW_;
    }
  }

  /**
  for (l=lmax; l>=1; l--){
    //    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )/sqrtK[l];
    phi_minus_1 = ( (2*l+1)*cotK*phi-phi_plus_1_times_sqrtK )*one_over_sqrtK[l];
    phi_plus_1_times_sqrtK = phi*sqrtK[l];
    phi = phi_minus_1;
    PhiL[l-1] = phi;
    if (fabs(phi)>_HYPER_OVERFLOW_){
      //
      phi *= _ONE_OVER_HYPER_OVERFLOW_;
      phi_plus_1_times_sqrtK *= _ONE_OVER_HYPER_OVERFLOW_;
      //Rescale whole Phi vector until this point:
      for (k=l-1; k<=lmax; k++)
        PhiL[k] *=_ONE_OVER_HYPER_OVERFLOW_;
    }
  }
  */
  scaling = phi0/phi;
  for (k=0; k<=lmax; k++)
    PhiL[k] *= scaling;
  return _SUCCESS_;
}

int hyperspherical_backwards_recurrence_chunk(int K,
                                              int lmax,
                                              double beta,
                                              double * __restrict__ x,
                                              double * __restrict__ sinK,
                                              double * __restrict__ cotK,
                                              int chunk,
                                              double * __restrict__ sqrtK,
                                              double * __restrict__ one_over_sqrtK,
                                              double * __restrict__ PhiL){
  double phi0, phi1, phipr1;
  int l, k, isign;
  int funcreturn = _FAILURE_;
  int index_x;
  double scalevec[_HYPER_CHUNK_]={0};

  for (index_x=0; index_x<chunk; index_x++){
    if (K==1){
      if (beta > 1.5*lmax) {
        funcreturn = get_CF1(K,lmax,beta,cotK[index_x], &phipr1, &isign);
      }
      if (funcreturn == _FAILURE_) {
        CF1_from_Gegenbauer(lmax,(int) (beta+0.2),sinK[index_x],cotK[index_x], &phipr1);
      }
      phi1 = 1.0;
    }
    else{
      get_CF1(K,lmax,beta,cotK[index_x], &phipr1, &isign);
      phi1 = isign;
      phipr1 *=phi1;
    }

    PhiL[lmax*chunk+index_x] = phi1;
    PhiL[(lmax-1)*chunk+index_x] = one_over_sqrtK[lmax]*
      ((lmax+1)*cotK[index_x]*phi1+phipr1);

  }
  for (l=lmax-2; l>=0; l--){
    //Use recurrence Phi_{l} = --Phi_{l+1} + -- Phi_{l+2}
    for (index_x=0; index_x<chunk; index_x++){
      PhiL[l*chunk+index_x] = one_over_sqrtK[l+1]*
        ((2*l+3)*cotK[index_x]*PhiL[(l+1)*chunk+index_x]-
         sqrtK[l+2]*PhiL[(l+2)*chunk+index_x]);
    }

    if (fabs(PhiL[l*chunk])>_HYPER_OVERFLOW_){
      //Rescale whole Phi vector until this point.
      //Create scale vector:
      for (index_x=0; index_x<chunk; index_x++)
        scalevec[index_x] = fabs(1.0/PhiL[l*chunk+index_x]);
      //Now do the scaling: (We do it this way to access elements in order)
      for (k=l; k<=lmax; k++){
        for (index_x=0; index_x<chunk; index_x++){
          PhiL[k*chunk+index_x] *= scalevec[index_x];
        }
      }
    }
  }
  for(index_x=0; index_x<chunk; index_x++){
    phi0 = sin(beta*x[index_x])/(beta*sinK[index_x]);
    scalevec[index_x] = phi0/PhiL[index_x];
  }
  for (k=0; k<=lmax; k++){
    for (index_x=0; index_x<chunk; index_x++){
      PhiL[k*chunk+index_x] *= scalevec[index_x];
    }
  }

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

  if (K==1) maxiter = (int)(beta-l-10);
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
  if (n<0)
    return _FAILURE_;
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
    G = 0.0;
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
  return _SUCCESS_;
}

 int hyperspherical_WKB_vec(int l,
                            double beta,
                            double * __restrict__ sinK_vec,
                            int size_sinK_vec,
                            double * __restrict__ Phi){
  double e, w, w2, alpha, alpha2, t;
  double S, Q, C, argu, Ai;
  int airy_sign = 1, phisign = 1;
  int index_sinK;
  double one_over_alpha;
  double one_over_alpha2;
  double one_over_sqrt_one_plus_alpha2;
  double sqrt_alpha;
  double one_over_e;
  double one_over_beta;
  double cscK;
  double pow_argu_onesixth;

  one_over_e = sqrt(l*(l+1.0));
  e = 1.0/one_over_e;
  alpha = beta*e;
  alpha2 = alpha*alpha;
  one_over_alpha = 1.0/alpha;
  one_over_alpha2 = one_over_alpha*one_over_alpha;
  one_over_sqrt_one_plus_alpha2 = 1.0/sqrt(1.0+alpha2);
  sqrt_alpha=sqrt(alpha);
  one_over_beta = 1.0/beta;

  for (index_sinK=0; index_sinK<size_sinK_vec; index_sinK++){
    cscK=1.0/sinK_vec[index_sinK];
    w = alpha*sinK_vec[index_sinK];
    w2 = w*w;
    if (alpha > cscK){
      S = alpha*log((sqrt(w2-1.0)+sqrt(w2+alpha2))*one_over_sqrt_one_plus_alpha2)+
        atan(one_over_alpha*sqrt((w2+alpha2)/(w2-1.0)))-M_PI_2;
      airy_sign = -1;
    }
    else{
      t = sqrt(1.0-w2)/sqrt(1.0+w2*one_over_alpha2);
      S = atanh(t)-alpha*atan(t*one_over_alpha);
      airy_sign = 1;
    }
    argu = 1.5*S*one_over_e;
    Q = cscK*cscK-alpha2;
    C = 0.5*sqrt_alpha*one_over_beta;
    pow_argu_onesixth = pow(argu,1.0/6.0);
    Ai = airy_cheb_approx(airy_sign*pow(pow_argu_onesixth,4));
    Phi[index_sinK] = phisign*2.0*_SQRT_PI_*C*pow_argu_onesixth*pow(fabs(Q),-0.25)*Ai*cscK;
  }
  return _SUCCESS_;
}


int hyperspherical_WKB(int K,int l,double beta,double y, double *Phi){
  double e, w, w2, alpha, alpha2, CscK, ytp, t;
  double S, Q, C, argu, Ai;
  int airy_sign = 1, phisign = 1, dphisign = 1, intbeta;
  double ldbl = l;

  if (K==1){
    //Limit range to [0; pi/2]:
    intbeta = (int)(beta+0.4); //Round to nearest integer (just to be sure)
    ClosedModY(l, intbeta, &y, &phisign, &dphisign);
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
    return _FAILURE_;
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
  else{
    return _FAILURE_;
  }
  argu = 3.0*S/(2.0*e);
  Q = CscK*CscK-alpha2;
  C = 0.5*sqrt(alpha)/beta;
  Ai = airy_cheb_approx(airy_sign*pow(argu,2.0/3.0));
  *Phi = phisign*2.0*sqrt(_PI_)*C*pow(argu,1.0/6.0)*pow(fabs(Q),-0.25)*Ai*CscK;
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
  const double A[5] = {1.1282427601,-0.6803534e-4,0.16687e-6,-0.128e-8,0.2e-10};
  const double B[5] = {0.7822108673e-1,-0.6895649e-4,0.32857e-6,-0.37e-8,0.7e-10};
  double x,y,t,Ai,zeta,theta,sintheta,costheta,FA,FB;

  x = -z;
  zeta = _TWO_OVER_THREE_*x*sqrt(x);
  theta = zeta+0.25*_PI_;
  sintheta = sin(theta);
  costheta = cos(theta);

  y = pow(7.0/x,3);//y = (7.0/x)*(7.0/x)*(7.0/x);
  FA = cheb(y,5,A);
  FB = cheb(y,5,B)/zeta;

  t = pow(x,-0.25);
  Ai = t*(sintheta*FA-costheta*FB);
  //Bi = t*(costheta*FA+sintheta*FB);
  return Ai;
}

double coef2(double z){
  const double A[17] = {0.11535880704,0.6542816649e-1,0.26091774326,0.21959346500,
              0.12476382168,-0.43476292594,0.28357718605,-0.9751797082e-1,
              0.2182551823e-1,-0.350454097e-2,0.42778312e-3,
              -0.4127564e-4,0.323880e-5,-0.21123e-6,0.1165e-7,
              -0.55e-9,0.2e-10};
  const double B[16] = {0.10888288487,-0.17511655051,0.13887871101,-0.11469998998,
             0.22377807641,-0.18546243714,0.8063565186e-1,
             -0.2208847864e-1,0.422444527e-2,-0.60131028e-3,
             0.6653890e-4,-0.590842e-5,0.43131e-6,-0.2638e-7,
             0.137e-8,-0.6e-10};
  //Ej = 3^(-j/3)/Gamma(j/3);
  double E1 = 0.355028053887817, E2 = 0.258819403792807;
  double x,FA,FB,Ai;
  x = -(z/7.0)*(z/7.0)*(z/7.0);
  FA = E1*cheb(x,17,A);
  FB = E2*z*cheb(x,16,B);
  Ai = FA-FB;
  //Bi = sqrt(3)*(FA+FB);
  return Ai;
}
double coef3(double z){
  const double A[20] = {1.2974695794,-0.20230907821,
              -0.45786169516,0.12953331987,0.6983827954e-1,
              -0.3005148746e-1,-0.493036334e-2,0.390425474e-2,
              -0.1546539e-4,-0.32193999e-3,0.3812734e-4,0.1714935e-4,
              -0.416096e-5,-0.50623e-6,0.26346e-6,-0.281e-8,
              -0.1122e-7,0.120e-8,0.31e-9,-0.7e-10};
  /**
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

  Ai = EY*cheb(x,20,A);
  //Bi = EX*cheb(x,25,&(B[0]));
  return Ai;
}

double coef4(double z){
  const double A[7]={0.56265126169,-0.76136219e-3,0.765252e-5,-0.14228e-6,
            0.380e-8,-0.13e-9,0.1e-10};
  /**  double B[7]={1.1316635302,0.166141673e-02,0.1968882e-04,0.47047e-06,
            0.1769e-7,0.94e-9,0.6e-10};
  */
  double x,Y,t,zeta,EX,EY,Ai;

  Y = z*sqrt(z);
  zeta = 2.0/3.0*Y;
  EX = exp(zeta);
  EY = 1.0/EX;
  x = 7*sqrt(7)/Y;
  t = pow(z,-0.25);

  Ai = t*EY*cheb(x, 7, A);
  //Bi = t*EX*cheb(x, 7, &(B[0]));
  return Ai;
}

double cheb(double x, int n, const double A[]){
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

int ClosedModY(int l, int beta, double *y, int * phisign, int * dphisign){

  *phisign = 1;
  *dphisign = 1;


  while (*y > _TWOPI_)
    *y -= _TWOPI_;

  if ((*y) > _PI_){
    *y = 2.0*_PI_-(*y);
    //phisign *= pow(-1,l)
    if (l%2==1) //l is odd
      *phisign = -*phisign;
    else        //l is even
      *dphisign = -*dphisign;
  }
  if ((*y)>0.5*_PI_){
    *y = _PI_-(*y);
    //phisign *= pow(-1,beta-l-1)
    if ((beta-l)%2==0) //beta-l-1 odd
      *phisign = -*phisign;
    else               //beta-l-1 even
      *dphisign = -*dphisign;
  }
  return _SUCCESS_;
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
                                                  index_l, x, Phi, NULL,NULL);
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

int hyperspherical_get_xmin_from_Airy(int K,
                                      int l,
                                      double beta,
                                      double xtol,
                                      double phiminabs,
                                      double *xmin,
                                      int *fevals){
  double xold, xtp=0, xleft, xright, xnew;
  double Fnew, Fold, Fleft, Fright;
  double delx, lambda;
  double AIRY_SAFETY = 1e-6;
  struct WKB_parameters wkbstruct;
  //Start searching from turning point:
  switch (K){
  case -1:
    xtp = asinh(sqrt(l*(l+1.))/beta);
    break;
  case 0:
    xtp = sqrt(l*(l+1.))/beta;
    break;
  case 1:
    xtp = asin(sqrt(l*(l+1.))/beta);
    break;
  }
  wkbstruct.K = K;
  wkbstruct.l = l;
  wkbstruct.beta = beta;
  wkbstruct.phiminabs = phiminabs;

  xnew = 0.99*xtp;

  Fnew = PhiWKB_minus_phiminabs(xnew,&wkbstruct);
  *fevals = (*fevals)+1;

  lambda = 2*_PI_/(beta+5.0);
  if (Fnew>0)
    delx = -lambda;
  else
    delx = 0.25*lambda;

  do {
    //printf("In the loop: xnew = %g, Fnew=%g, Fold=%g\n",xnew,Fnew,Fold);
    xold = xnew;
    Fold = Fnew;
    xnew += delx;
    if (xnew<AIRY_SAFETY){
      xnew = AIRY_SAFETY;
      Fnew = PhiWKB_minus_phiminabs(xnew,&wkbstruct);
      *fevals = (*fevals)+1;
      if (Fnew>=0.0){
        *xmin = xnew;
        return _SUCCESS_;
      }
      else{
        break;
      }
    }
    Fnew = PhiWKB_minus_phiminabs(xnew,&wkbstruct);
    *fevals = (*fevals)+1;
  } while (SIGN(Fnew)==(SIGN(Fold)));

  if (Fnew<=0.0){
    xleft = xnew;
    Fleft = Fnew;
    xright = xold;
    Fright = Fold;
  }
  else{
    xleft = xold;
    Fleft = Fold;
    xright = xnew;
    Fright = Fnew;
  }

  fzero_ridder(PhiWKB_minus_phiminabs,
               xleft,
               xright,
               xtol,
               &wkbstruct,
               &Fleft,
               &Fright,
               xmin,
               fevals);

  return _SUCCESS_;
}

double PhiWKB_minus_phiminabs(double x, void *param){
   double phiwkb;
   struct WKB_parameters *wkbparam = (struct WKB_parameters *)param;
   hyperspherical_WKB(wkbparam->K,wkbparam->l,wkbparam->beta,x, &phiwkb);
   return(fabs(phiwkb)-wkbparam->phiminabs);
}

int fzero_ridder(double (*func)(double, void *),
                  double x1,
                  double x2,
                  double xtol,
                  void *param,
                  double *Fx1,
                  double *Fx2,
                  double *xzero,
                  int *fevals){
   /**Using Ridders' method, return the root of a function func known to
      lie between x1 and x2. The root, returned as zriddr, will be found to
      an approximate accuracy xtol.
   */
   int j,MAXIT=1000;
    double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
    if ((Fx1!=NULL)&&(Fx2!=NULL)){
      fl = *Fx1;
      fh = *Fx2;
    }
    else{
      fl=(*func)(x1,param);
      fh=(*func)(x2,param);
      *fevals = (*fevals)+2;
    }
    if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
      xl=x1;
      xh=x2;
      ans=-1.11e11;
      for (j=1;j<=MAXIT;j++) {
        xm=0.5*(xl+xh);
        fm=(*func)(xm,param);
        *fevals = (*fevals)+1;
        s=sqrt(fm*fm-fl*fh);
        if (s == 0.0){
          *xzero = ans;
          //printf("Success 1\n");
          return _SUCCESS_;
        }
        xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
        if (fabs(xnew-ans) <= xtol) {
          *xzero = ans;
          return _SUCCESS_;
        }
        ans=xnew;
        fnew=(*func)(ans,param);
        *fevals = (*fevals)+1;
        if (fnew == 0.0){
          *xzero = ans;
          //printf("Success 2, ans=%g\n",ans);
          return _SUCCESS_;
        }
        if (NRSIGN(fm,fnew) != fm) {
          xl=xm;
          fl=fm;
          xh=ans;
          fh=fnew;
        } else if (NRSIGN(fl,fnew) != fl) {
          xh=ans;
          fh=fnew;
        } else if (NRSIGN(fh,fnew) != fh) {
          xl=ans;
          fl=fnew;
        } else return _FAILURE_;
        if (fabs(xh-xl) <= xtol) {
          *xzero = ans;
          //        printf("Success 3\n");
          return _SUCCESS_;
        }
      }
      fprintf(stderr,"zriddr exceed maximum iterations");
      return _FAILURE_;
    }
    else {
      if (fl == 0.0) return x1;
      if (fh == 0.0) return x2;
      fprintf(stderr,"root must be bracketed in zriddr.");
      return _FAILURE_;
    }
    fprintf(stderr,"Failure in Ridder\n");
    return _FAILURE_;
  }

int HypersphericalExplicit(int K,int l, double beta,double x, double *Phi){
   /** Explicit formulae for the Hyperspherical Besselfunctions of order
l<=9.
       phi_tilde = gam * beta * cos(x*beta) + delta * sin(x*beta),
       and Phi = phi_tilde *cscK/sqrt(NK). Gamma and delta are
polynomials in
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
     NK = beta2*(beta2*(202759531776.0 + 32946*beta12 + beta16 + 15088541896.0*beta4 + 68943381.0*beta8) -
         5*(26336378880.0 + 19*beta4*(893321712.0 + 3*beta12 + 14395719.0*beta4 + 21046.0*beta8))*K);
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

int hyperspherical_get_xmin_from_approx(int K,
                                        int l,
                                        double nu,
                                        double ignore1,
                                        double phiminabs,
                                        double *xmin,
                                        int *ignore2){

  double l_plus_half;
  double lhs;
  double alpha;
  double ldbl = l;
  double x;

  l_plus_half = l+0.5;
  lhs = 1.0/l_plus_half*log(2*phiminabs*l_plus_half);
  //Using Chebyshev cubic root, much cleaner:
  alpha = -2.0*lhs/5.0*(1.0+2.0*cosh(1.0/3.0*acosh(1.0+375.0/(16.0*lhs*lhs))));
  x = l_plus_half/cosh(alpha)/nu;
  if (K==-1){
    //%Correct for open case:
    x *= asinh(ldbl/nu)/(ldbl/nu);
    //...and fudge for small nu:
    x *=((nu+0.4567)/(nu+1.24)-2.209e-3);
  }
  else if(K==1){
    //Correct for closed case if possible
    x *= asin(ldbl/nu)/(ldbl/nu);
  }
  *xmin = x;
  return _SUCCESS_;
}

/** Generate the 2^3-1 non-trivial versions of the functions
    hyperspherical_Hermite3_interpolation_vectorXXX(),
    hyperspherical_Hermite4_interpolation_vectorXXX() and
    hyperspherical_Hermite6_interpolation_vectorXXX() using the
    preprocessor. Apologise in advance, but speed for this function
    is important and it is better than manual copy-paste.
*/
int hyperspherical_Hermite3_interpolation_vector_Phi(HyperInterpStruct *pHIS,
                                                     int nxi,
                                                     int lnum,
                                                     double * xinterp,
                                                     double * Phi,
                                                     ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite3_interpolation_vector_dPhi(HyperInterpStruct *pHIS,
                                                      int nxi,
                                                      int lnum,
                                                      double * xinterp,
                                                      double * dPhi,
                                                      ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_DPHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite3_interpolation_vector_d2Phi(HyperInterpStruct *pHIS,
                                                       int nxi,
                                                       int lnum,
                                                       double * xinterp,
                                                       double * d2Phi,
                                                       ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_D2PHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite3_interpolation_vector_PhidPhi(HyperInterpStruct *pHIS,
                                                         int nxi,
                                                         int lnum,
                                                         double * xinterp,
                                                         double * Phi,
                                                         double * dPhi,
                                                         ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_DPHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite3_interpolation_vector_Phid2Phi(HyperInterpStruct *pHIS,
                                                          int nxi,
                                                          int lnum,
                                                          double * xinterp,
                                                          double * Phi,
                                                          double * d2Phi,
                                                          ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_D2PHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite3_interpolation_vector_dPhid2Phi(HyperInterpStruct *pHIS,
                                                           int nxi,
                                                           int lnum,
                                                           double * xinterp,
                                                           double * dPhi,
                                                           double * d2Phi,
                                                           ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_DPHI
#define HERMITE_DO_D2PHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}

int hyperspherical_Hermite3_interpolation_vector_PhidPhid2Phi(HyperInterpStruct *pHIS,
                                                              int nxi,
                                                              int lnum,
                                                              double * xinterp,
                                                              double *Phi,
                                                              double * dPhi,
                                                              double * d2Phi,
                                                              ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_DPHI
#define HERMITE_DO_D2PHI
#include "hermite3_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite4_interpolation_vector_Phi(HyperInterpStruct *pHIS,
                                                     int nxi,
                                                     int lnum,
                                                     double * xinterp,
                                                     double * Phi,
                                                     ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite4_interpolation_vector_dPhi(HyperInterpStruct *pHIS,
                                                      int nxi,
                                                      int lnum,
                                                      double * xinterp,
                                                      double * dPhi,
                                                      ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_DPHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite4_interpolation_vector_d2Phi(HyperInterpStruct *pHIS,
                                                       int nxi,
                                                       int lnum,
                                                       double * xinterp,
                                                       double * d2Phi,
                                                       ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_D2PHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite4_interpolation_vector_PhidPhi(HyperInterpStruct *pHIS,
                                                         int nxi,
                                                         int lnum,
                                                         double * xinterp,
                                                         double * Phi,
                                                         double * dPhi,
                                                         ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_DPHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite4_interpolation_vector_Phid2Phi(HyperInterpStruct *pHIS,
                                                          int nxi,
                                                          int lnum,
                                                          double * xinterp,
                                                          double * Phi,
                                                          double * d2Phi,
                                                          ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_D2PHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite4_interpolation_vector_dPhid2Phi(HyperInterpStruct *pHIS,
                                                           int nxi,
                                                           int lnum,
                                                           double * xinterp,
                                                           double * dPhi,
                                                           double * d2Phi,
                                                           ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_DPHI
#define HERMITE_DO_D2PHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}

int hyperspherical_Hermite4_interpolation_vector_PhidPhid2Phi(HyperInterpStruct *pHIS,
                                                              int nxi,
                                                              int lnum,
                                                              double * xinterp,
                                                              double *Phi,
                                                              double * dPhi,
                                                              double * d2Phi,
                                                              ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_DPHI
#define HERMITE_DO_D2PHI
#include "hermite4_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite6_interpolation_vector_Phi(HyperInterpStruct *pHIS,
                                                     int nxi,
                                                     int lnum,
                                                     double * xinterp,
                                                     double * Phi,
                                                     ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite6_interpolation_vector_dPhi(HyperInterpStruct *pHIS,
                                                      int nxi,
                                                      int lnum,
                                                      double * xinterp,
                                                      double * dPhi,
                                                      ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_DPHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite6_interpolation_vector_d2Phi(HyperInterpStruct *pHIS,
                                                       int nxi,
                                                       int lnum,
                                                       double * xinterp,
                                                       double * d2Phi,
                                                       ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_D2PHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite6_interpolation_vector_PhidPhi(HyperInterpStruct *pHIS,
                                                         int nxi,
                                                         int lnum,
                                                         double * xinterp,
                                                         double * Phi,
                                                         double * dPhi,
                                                         ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_DPHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite6_interpolation_vector_Phid2Phi(HyperInterpStruct *pHIS,
                                                          int nxi,
                                                          int lnum,
                                                          double * xinterp,
                                                          double * Phi,
                                                          double * d2Phi,
                                                          ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_D2PHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}
int hyperspherical_Hermite6_interpolation_vector_dPhid2Phi(HyperInterpStruct *pHIS,
                                                           int nxi,
                                                           int lnum,
                                                           double * xinterp,
                                                           double * dPhi,
                                                           double * d2Phi,
                                                           ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_DPHI
#define HERMITE_DO_D2PHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}

int hyperspherical_Hermite6_interpolation_vector_PhidPhid2Phi(HyperInterpStruct *pHIS,
                                                              int nxi,
                                                              int lnum,
                                                              double * xinterp,
                                                              double *Phi,
                                                              double * dPhi,
                                                              double * d2Phi,
                                                              ErrorMsg error_message) {
#undef HERMITE_DO_PHI
#undef HERMITE_DO_DPHI
#undef HERMITE_DO_D2PHI
#define HERMITE_DO_PHI
#define HERMITE_DO_DPHI
#define HERMITE_DO_D2PHI
#include "hermite6_interpolation_csource.h"
  return _SUCCESS_;
}

