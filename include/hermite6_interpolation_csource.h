/** Hermite interpolation of order 6 for Phi, dPhi, and d2Phi. When xinterp
    is sorted (increasing), computations can be reused. On the other hand,
    for a randomly called value, the routine is not much slower than a
    routine optimised for this case. The more sorted the vector, the faster
    the execution time. For closed case, the interpolation structure only
    covers [safety;pi/2-safety]. The calling routine should respect this.
    if sinK and cosK are not NULL, we will also interpolate them.
*/

double ym=0, yp=0, dym=0, dyp=0, d2ym=0, d2yp=0, x, z, z2, z3, z4, z5;
double cotKm=0,cotKp=0,sinKm=0,sinKp=0, sinKm2, sinKp2;
#ifdef HERMITE_DO_PHI
double a1=0, a2=0, a3=0, a4=0, a5=0;
#endif
#ifdef HERMITE_DO_DPHI
double b1=0, b2=0, b3=0, b4=0, b5=0;
#endif
#ifdef HERMITE_DO_D2PHI
double c1=0, c2=0, c3=0, c4=0, c5=0;
double d4ym=0, d4yp=0;
#endif
#if defined (HERMITE_DO_DPHI) || defined (HERMITE_DO_D2PHI)
double d3ym = 0, d3yp=0;
#endif
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
#ifdef HERMITE_DO_PHI
    Phi[j] = 0.0;
#endif
#ifdef HERMITE_DO_DPHI
    dPhi[j] = 0.0;
#endif
#ifdef HERMITE_DO_D2PHI
    d2Phi[j] = 0.0;
#endif
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
#if defined HERMITE_DO_DPHI || defined HERMITE_DO_D2PHI
      d3ym = -2*cotKm*d2ym-2*ym*lxlp1*cotKm/sinKm2+
        dym*(K-beta2+(2+lxlp1)/sinKm2);
#endif
#ifdef HERMITE_DO_D2PHI
      d4ym = -2*cotKm*d3ym + d2ym*(K-beta2+(4+lxlp1)/sinKm2)+
        dym*(-4*(1+lxlp1)*cotKm/sinKm2)+
        ym*(2*lxlp1/sinKm2*(2*cotKm*cotKm+1/sinKm2));
#endif
    }
    else{
      //x>current_border but not next border: I have moved to next block.
      current_border_idx++;
      //printf("Current border index at else: %d\n",current_border_idx);
      //Copy former right derivatives to left derivatives.
      ym = yp;
      dym = dyp;
      d2ym = d2yp;
#if defined HERMITE_DO_DPHI || defined HERMITE_DO_D2PHI
      d3ym = d3yp;
#endif
#ifdef HERMITE_DO_D2PHI
      d4ym = d4yp;
#endif
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
#if defined HERMITE_DO_DPHI || defined HERMITE_DO_D2PHI
    d3yp = -2*cotKp*d2yp-2*yp*lxlp1*cotKp/sinKp2+
      dyp*(K-beta2+(2+lxlp1)/sinKp2);
#endif
#ifdef HERMITE_DO_D2PHI
    d4yp = -2*cotKp*d3yp + d2yp*(K-beta2+(4+lxlp1)/sinKp2)+
      dyp*(-4*(1+lxlp1)*cotKp/sinKp2)+
      yp*(2*lxlp1/sinKp2*(2*cotKp*cotKp+1/sinKp2));
#endif

#ifdef HERMITE_DO_PHI
    a1 = dym*deltax;
    a2 = 0.5*d2ym*deltax2;
    a3 = (-1.5*d2ym+0.5*d2yp)*deltax2-(6*dym+4*dyp)*deltax-10*(ym-yp);
    a4 = (1.5*d2ym-d2yp)*deltax2+(8*dym+7*dyp)*deltax+15*(ym-yp);
    a5 = (-0.5*d2ym+0.5*d2yp)*deltax2-3*(dym+dyp)*deltax-6*(ym-yp);
#endif
#ifdef HERMITE_DO_DPHI
    b1 = d2ym*deltax;
    b2 = 0.5*d3ym*deltax2;
    b3 = (-1.5*d3ym+0.5*d3yp)*deltax2-(6*d2ym+4*d2yp)*deltax-10*(dym-dyp);
    b4 = (1.5*d3ym-d3yp)*deltax2+(8*d2ym+7*d2yp)*deltax+15*(dym-dyp);
    b5 = (-0.5*d3ym+0.5*d3yp)*deltax2-3*(d2ym+d2yp)*deltax-6*(dym-dyp);
#endif
#ifdef HERMITE_DO_D2PHI
    c1 = d3ym*deltax;
    c2 = 0.5*d4ym*deltax2;
    c3 = (-1.5*d4ym+0.5*d4yp)*deltax2-(6*d3ym+4*d3yp)*deltax-10*(d2ym-d2yp);
    c4 = (1.5*d4ym-d4yp)*deltax2+(8*d3ym+7*d3yp)*deltax+15*(d2ym-d2yp);
    c5 = (-0.5*d4ym+0.5*d4yp)*deltax2-3*(d3ym+d3yp)*deltax-6*(d2ym-d2yp);
#endif
  }
  //Evaluate polynomial:
  z = (x-left_border)/deltax;
  z2 = z*z;
  z3 = z2*z;
  z4 = z2*z2;
  z5 = z2*z3;
#ifdef HERMITE_DO_PHI
    Phi[j] = (ym+a1*z+a2*z2+a3*z3+a4*z4+a5*z5)*phisign;
#endif
#ifdef HERMITE_DO_DPHI
    dPhi[j] = (dym+b1*z+b2*z2+b3*z3+b4*z4+b5*z5)*dphisign;
#endif
#ifdef HERMITE_DO_D2PHI
    d2Phi[j] = (d2ym+c1*z+c2*z2+c3*z3+c4*z4+c5*z5)*phisign;
#endif
 }

