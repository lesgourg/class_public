/** Hermite interpolation of order 4 for Phi, dPhi, and d2Phi. When xinterp
    is sorted (increasing), computations can be reused. On the other hand,
    for a randomly called value, the routine is not much slower than a
    routine optimised for this case. The more sorted the vector, the faster
    the execution time. For closed case, the interpolation structure only
    covers [safety;pi/2-safety]. The calling routine should respect this.
    if sinK and cosK are not NULL, we will also interpolate them.
*/

int l = pHIS->l[lnum];
double ym=0, yp=0, dym=0, dyp=0, x;
double z[3]={0.,0.,0.};
#ifdef HERMITE_DO_PHI
double a[3]={0.,0.,0.};
#endif
#ifdef HERMITE_DO_DPHI
double b[3]={0.,0.,0.};
#endif
#ifdef HERMITE_DO_D2PHI
double c[3]={0.,0.,0.};
double d3ym=0, d3yp=0;
#endif
#if defined (HERMITE_DO_DPHI) || defined (HERMITE_DO_D2PHI)
double *sinK = pHIS->sinK;
double *cotK = pHIS->cotK;
double cotKm=0,cotKp=0,sinKm=0,sinKp=0;
double sinKm2, sinKp2;
double d2ym = 0, d2yp=0;
int K = pHIS->K;
double lxlp1 = l*(l+1.0);
double beta = pHIS->beta;
double beta2 = beta*beta;
#endif
double *xvec;
double xmin, xmax, deltax;
double left_border, right_border, next_border;
int j, nx, current_border_idx=0;
double *Phi_l, *dPhi_l;
int phisign = 1, dphisign = 1;

/** Set logical flags. The compiler should probably generate 2^3-1=7
    different functions, according to these flags. If not, maybe I should
    do it.
*/

xvec = pHIS->x;
deltax = pHIS->delta_x;
nx = pHIS->x_size;
Phi_l = pHIS->phi+lnum*nx;
dPhi_l = pHIS->dphi+lnum*nx;

xmin = xvec[0];
xmax = xvec[nx-1];

left_border = xmax;
right_border = xmin;
next_border = xmin;

for (j=0; j<nxi; j++){
  x = xinterp[j];
  //take advantage of periodicity of functions in closed case
  if (pHIS->K==1)
    ClosedModY(l, (int)(pHIS->beta+0.2), &x, &phisign, &dphisign);
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
      ym = Phi_l[current_border_idx-1];
      dym = dPhi_l[current_border_idx-1];
#if defined HERMITE_DO_DPHI || defined HERMITE_DO_D2PHI
      cotKm = cotK[current_border_idx-1];
      sinKm = sinK[current_border_idx-1];
      sinKm2 = sinKm*sinKm;
      d2ym = -2*dym*cotKm+ym*(lxlp1/sinKm2-beta2+K);
#endif
#ifdef HERMITE_DO_D2PHI
      d3ym = -2*cotKm*d2ym-2*ym*lxlp1*cotKm/sinKm2+
        dym*(K-beta2+(2+lxlp1)/sinKm2);
#endif
    }
    else{
      //x>current_border but not next border: I have moved to next block.
      current_border_idx++;
      //printf("Current border index at else: %d\n",current_border_idx);
      //Copy former right derivatives to left derivatives.
      ym = yp;
      dym = dyp;
#if defined HERMITE_DO_DPHI || defined HERMITE_DO_D2PHI
      d2ym = d2yp;
      sinKm = sinKp;
      cotKm = cotKp;
#endif
#ifdef HERMITE_DO_D2PHI
      d3ym = d3yp;
#endif
    }
    left_border = xvec[MAX(0,current_border_idx-1)];
    right_border = xvec[current_border_idx];
    next_border = xvec[MIN(nx-1,current_border_idx+1)];
    //Evaluate right derivatives and calculate coefficients:
    yp = Phi_l[current_border_idx];
    dyp = dPhi_l[current_border_idx];
#if defined HERMITE_DO_DPHI || defined HERMITE_DO_D2PHI
    cotKp = cotK[current_border_idx];
    sinKp = sinK[current_border_idx];
    sinKp2 = sinKp*sinKp;
    d2yp = -2*dyp*cotKp+yp*(lxlp1/sinKp2-beta2+K);
#endif
#ifdef HERMITE_DO_D2PHI
    d3yp = -2*cotKp*d2yp-2*yp*lxlp1*cotKp/sinKp2+
      dyp*(K-beta2+(2+lxlp1)/sinKp2);
#endif

#ifdef HERMITE_DO_PHI
    a[0] = dym*deltax;
    a[1] = -2*dym*deltax-dyp*deltax-3*ym+3*yp;
    a[2] = dym*deltax+dyp*deltax+2*ym-2*yp;
#endif
#ifdef HERMITE_DO_DPHI
    b[0] = d2ym*deltax;
    b[1] = -2*d2ym*deltax-d2yp*deltax-3*dym+3*dyp;
    b[2] = d2ym*deltax+d2yp*deltax+2*dym-2*dyp;
#endif
#ifdef HERMITE_DO_D2PHI
    c[0] = d3ym*deltax;
    c[1] = -2*d3ym*deltax-d3yp*deltax-3*d2ym+3*d2yp;
    c[2] = d3ym*deltax+d3yp*deltax+2*d2ym-2*d2yp;
#endif
  }
  //Evaluate polynomial:
  z[0] = (x-left_border)/deltax;
  z[1] = z[0]*z[0];
  z[2] = z[1]*z[0];
#ifdef HERMITE_DO_PHI
  Phi[j] = (ym+a[0]*z[0]+a[1]*z[1]+a[2]*z[2])*phisign;
#endif
#ifdef HERMITE_DO_DPHI
  dPhi[j] = (dym+b[0]*z[0]+b[1]*z[1]+b[2]*z[2])*dphisign;
#endif
#ifdef HERMITE_DO_D2PHI
  d2Phi[j] = (d2ym+c[0]*z[0]+c[1]*z[1]+c[2]*z[2])*phisign;
#endif
 }

