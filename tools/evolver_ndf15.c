/****************************************/
/* Stiff ODE evolver for CLASS                    */
/* 19/11 2010                                     */
/* Thomas Tram                                    */
/****************************************/
/*    This is an variable order, adaptive stepsize ODE evolver for CLASS.
    It is based on the Numerical Differentiaion Formulas    of order 1-5,
    and can be used to solve stiff problems. The algorithm is described in
    [The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt, SIAM Journal
    on Scientific Computing, 18-1, 1997].

    The code will call the (*output) routine at x-values in t_vec[]. It
    will interpolate to get the y-value aswell as (optionally) the first and
    second derivative at points returned by the routine relevant_indices.
    Since the transfer function only depends on a few of these values, there
    is no need to interpolate all of them.

    Every time there is a stepsize change, we need to rebuild **dif, the matrix
    which holds the backward differences, to reflect the new stepsize. This is done
    in "adjust_stepsize" and is essentially a matrix multiplication. Every times
    there is either a stepsize change, the method order k is changed, or we compute
    a new Jacobian, we must call "new_linearisation" to calculate a new LU
    decomposition of a matrix.

    The method will not recompute the Jacobian in every step, but only when the
    Newton iterations fail to converge fast enough. This feature makes the
    solver competitive even for non-stiff problems.

    Statistics is saved in the stepstat[6] vector. The entries are:
    stepstat[0] = Successful steps.
    stepstat[1] = Failed steps.
    stepstat[2] = Total number of function evaluations.
    stepstat[3] = Number of Jacobians computed.
    stepstat[4] = Number of LU decompositions.
    stepstat[5] = Number of linear solves.
    If ppt->perturbations_verbose > 2, this statistic is printed at the end of
    each call to evolver.

    Sparsity:
    When the number of equations becomes high, too much times is spent on solving
    linear equations. Since the ODE-functions are not very coupled, one would normally
    supply a sparsity pattern of the jacobian, which would encode the couplings of the
    equations. However, this would be of some inconvenience to the users who just want
    to add some new physics without thinking too much about how the code work. So we
    pay a few more function evaluations, and calculate the full jacobian every time.

    Then, if jac->use_sparse==_TRUE_, numjac will try to construct a sparse matrix from
    the dense matrix. If there are too many nonzero elements in the dense matrix, numjac
    will stop constructing the sparse matrix and set jac->use_sparse=_FALSE_. The sparse
    matrix is stored in the compressed column format. (See sparse.h).

    In the sparse case, we also do partial pivoting, but with diagonal preference. The
    structure of the equations are nearly optimal for the LU decomposition, so we don't
    want to mess it up by too many row permutations if we can avoid it. This is also why
    do not use any column permutation to pre-order the matrix.
*/
#include "common.h"
#include "evolver_ndf15.h"
//#include "perturbations.h"
#include "sparse.h"

int evolver_ndf15(
          int (*derivs)(double x,double * y,double * dy,
                void * parameters_and_workspace, ErrorMsg error_message),
          double x_ini,
          double x_final,
          double * y_inout,
          int * used_in_output,
          int neq,
          void * parameters_and_workspace_for_derivs,
          double rtol,
          double minimum_variation,
          int (*timescale_and_approximation)(double x,
                             void * parameters_and_workspace,
                             double * timescales,
                             ErrorMsg error_message),
          double timestep_over_timescale,
          double * t_vec,
          int tres,
          int (*output)(double x,double y[],double dy[],int index_x,void * parameters_and_workspace,
                ErrorMsg error_message),
          int (*print_variables)(double x, double y[], double dy[], void *parameters_and_workspace,
                     ErrorMsg error_message),
          ErrorMsg error_message){

  /* Constants: */
  double G[5]={1.0,3.0/2.0,11.0/6.0,25.0/12.0,137.0/60.0};
  double alpha[5]={-37.0/200,-1.0/9.0,-8.23e-2,-4.15e-2, 0};
  double invGa[5],erconst[5];
  double abstol = 1e-15, eps=1e-16, threshold=abstol;
  int maxit=4, maxk=5;

  /* Logicals: */
  int Jcurrent,havrate,done,at_hmin,nofailed,gotynew,tooslow,*interpidx;

  /* Storage: */
  double *f0,*y,*wt,*ddfddt,*pred,*ynew,*invwt,*rhs,*psi,*difkp1,*del,*yinterp;
  double *tempvec1,*tempvec2,*ypinterp,*yppinterp;
  double **dif;
  struct jacobian jac;
  struct numjac_workspace nj_ws;

  /* Method variables: */
  double t,t0,tfinal,tnew=0;
  double rh,htspan,absh,hmin,hmax,h,tdel;
  double abshlast,hinvGak,minnrm,oldnrm=0.,newnrm;
  double err,hopt,errkm1,hkm1,errit,rate=0.,temp,errkp1,hkp1,maxtmp;
  int k,klast,nconhk,iter,next,kopt,tdir;

  /* Misc: */
  int stepstat[6],nfenj,j,ii,jj, numidx, neqp=neq+1;
  int verbose=0;
  int funcreturn;

  /** Allocate memory . */

  void * buffer;
  int buffer_size;

  buffer_size=
    15*neqp*sizeof(double)
    +neqp*sizeof(int)
    +neqp*sizeof(double*)
    +(7*neq+1)*sizeof(double);

  class_alloc(buffer,
          buffer_size,
          error_message);

  f0       =(double*)buffer;
  wt       =f0+neqp;
  ddfddt   =wt+neqp;
  pred     =ddfddt+neqp;
  y        =pred+neqp;
  invwt    =y+neqp;
  rhs      =invwt+neqp;
  psi      =rhs+neqp;
  difkp1   =psi+neqp;
  del      =difkp1+neqp;
  yinterp  =del+neqp;
  ypinterp =yinterp+neqp;
  yppinterp=ypinterp+neqp;
  tempvec1 =yppinterp+neqp;
  tempvec2 =tempvec1+neqp;

  interpidx=(int*)(tempvec2+neqp);

  dif      =(double**)(interpidx+neqp);
  dif[1]   =(double*)(dif+neqp);
  for(j=2;j<=neq;j++) dif[j] = dif[j-1]+7; /* Set row pointers... */
  dif[0] = NULL;
  /* for (ii=0;ii<(7*neq+1);ii++) dif[1][ii]=0.; */
  for (j=1; j<=neq; j++) {
    for (ii=1;ii<=7;ii++) {
      dif[j][ii]=0.;
    }
  }

  /*     class_alloc(f0,sizeof(double)*neqp,error_message); */
  /*     class_alloc(wt,sizeof(double)*neqp,error_message); */
  /*     class_alloc(ddfddt,sizeof(double)*neqp,error_message); */
  /*     class_alloc(pred,sizeof(double)*neqp,error_message); */
  /*     class_alloc(y,sizeof(double)*neqp,error_message); */
  /*     class_alloc(invwt,sizeof(double)*neqp,error_message); */
  /*     class_alloc(rhs,sizeof(double)*neqp,error_message); */
  /*     class_alloc(psi,sizeof(double)*neqp,error_message); */
  /*     class_alloc(difkp1,sizeof(double)*neqp,error_message); */
  /*     class_alloc(del,sizeof(double)*neqp,error_message); */
  /*     class_alloc(yinterp,sizeof(double)*neqp,error_message); */
  /*     class_alloc(ypinterp,sizeof(double)*neqp,error_message); */
  /*     class_alloc(yppinterp,sizeof(double)*neqp,error_message); */
  /*     class_alloc(tempvec1,sizeof(double)*neqp,error_message); */
  /*     class_alloc(tempvec2,sizeof(double)*neqp,error_message); */

  /*     class_alloc(interpidx,sizeof(int)*neqp,error_message); */

  /* Allocate vector of pointers to rows of dif:*/
  /*     class_alloc(dif,sizeof(double*)*neqp,error_message);  */
  /*     class_calloc(dif[1],(7*neq+1),sizeof(double),error_message); */
  /*     dif[0] = NULL; */
  /*     for(j=2;j<=neq;j++) dif[j] = dif[j-1]+7; */ /* Set row pointers... */

  /*Set pointers:*/
  ynew = y_inout-1; /* This way y_inout is always up to date. */

  /*Initialize the jacobian:*/
  class_call(initialize_jacobian(&jac,neq,error_message),error_message,error_message);

  /* Initialize workspace for numjac: */
  class_call(initialize_numjac_workspace(&nj_ws,neq,error_message),error_message,error_message);

  /* Initialize some method parameters:*/
  for(ii=0;ii<5;ii++){
    invGa[ii] = 1.0/(G[ii]*(1.0 - alpha[ii]));
    erconst[ii] = alpha[ii]*G[ii] + 1.0/(2.0+ii);
  }

  /* Set the relevant indices which needs to be found by interpolation. */
  /* But if we want to print variables for testing purposes, just interpolate everything.. */
  for(ii=1;ii<=neq;ii++){
    y[ii] = y_inout[ii-1];
    if (print_variables==NULL){
      interpidx[ii]=used_in_output[ii-1];
    }
    else{
      interpidx[ii]=1;
    }
  }
  t0 = x_ini;
  tfinal = x_final;

  /* Some CLASS-specific stuff:*/
  next=0;
  while (t_vec[next] < t0) next++;

  if (verbose > 1){
    numidx=0;
    for(ii=1;ii<=neq;ii++){
      if (interpidx[ii]==_TRUE_) numidx++;
    }
    printf("%d/%d ",numidx,neq);
  }

  htspan = fabs(tfinal-t0);
  for(ii=0;ii<6;ii++) stepstat[ii] = 0;

  class_call((*derivs)(t0,y+1,f0+1,parameters_and_workspace_for_derivs,error_message),error_message,error_message);
  stepstat[2] +=1;
  if ((tfinal-t0)<0.0){
    tdir = -1;
  }
  else{
    tdir = 1;
  }
  hmax = (tfinal-t0)/10.0;
  t = t0;


  nfenj=0;
  class_call(numjac((*derivs),t,y,f0,&jac,&nj_ws,abstol,neq,
             &nfenj,parameters_and_workspace_for_derivs,error_message),
             error_message,error_message);
  stepstat[3] += 1;
  stepstat[2] += nfenj;
  Jcurrent = _TRUE_; /* True */

  hmin = 16.0*eps*MAX(fabs(t),fabs(tfinal));
  /*Calculate initial step */
  rh = 0.0;

  for(jj=1;jj<=neq;jj++){
    wt[jj] = MAX(fabs(y[jj]),threshold);
    /*printf("wt: %4.8f \n",wt[jj]);*/
    rh = MAX(rh,1.25/sqrt(rtol)*fabs(f0[jj]/wt[jj]));
  }

  absh = MIN(hmax, htspan);
  if (absh * rh > 1.0) absh = 1.0 / rh;

  absh = MAX(absh, hmin);
  h = tdir * absh;
  tdel = (t + tdir*MIN(sqrt(eps)*MAX(fabs(t),fabs(t+h)),absh)) - t;

  class_call((*derivs)(t+tdel,y+1,tempvec1+1,parameters_and_workspace_for_derivs,error_message),
             error_message,error_message);
  stepstat[2] += 1;

  /*I assume that a full jacobi matrix is always calculated in the beginning...*/
  for(ii=1;ii<=neq;ii++){
    ddfddt[ii]=0.0;
    for(jj=1;jj<=neq;jj++){
      ddfddt[ii]+=(jac.dfdy[ii][jj])*f0[jj];
    }
  }

  rh = 0.0;
  for(ii=1;ii<=neq;ii++){
    ddfddt[ii] += (tempvec1[ii] - f0[ii]) / tdel;
    rh = MAX(rh,1.25*sqrt(0.5*fabs(ddfddt[ii]/wt[ii])/rtol));
  }
  absh = MIN(hmax, htspan);
  if (absh * rh > 1.0) absh = 1.0 / rh;

  absh = MAX(absh, hmin);
  h = tdir * absh;
  /* Done calculating initial step
     Get ready to do the loop:*/
  k = 1;            /*start at order 1 with BDF1    */
  klast = k;
  abshlast = absh;

  for(ii=1;ii<=neq;ii++) dif[ii][1] = h*f0[ii];

  hinvGak = h*invGa[k-1];
  nconhk = 0;     /*steps taken with current h and k*/
  class_call(new_linearisation(&jac,hinvGak,neq,error_message),
             error_message,error_message);
  stepstat[4] += 1;
  havrate = _FALSE_; /*false*/

  /* Doing main loop: */
  done = _FALSE_;
  at_hmin = _FALSE_;
  while (done==_FALSE_){
    /**class_test(stepstat[2] > 1e5, error_message,
           "Too many steps in evolver! Current stepsize:%g, in interval: [%g:%g]\n",
           absh,t0,tfinal);*/
    maxtmp = MAX(hmin,absh);
    absh = MIN(hmax, maxtmp);
    if (fabs(absh-hmin)<100*eps){
      /* If the stepsize has not changed */
      if (at_hmin==_TRUE_){
    absh = abshlast;    /*required by stepsize recovery */
      }
      at_hmin = _TRUE_;
    }
    else{
      at_hmin = _FALSE_;
    }
    h = tdir * absh;
    /* Stretch the step if within 10% of tfinal-t. */
    if (1.1*absh >= fabs(tfinal - t)){
      h = tfinal - t;
      absh = fabs(h);
      done = _TRUE_;
    }
    if (((fabs(absh-abshlast)/absh)>1e-6)||(k!=klast)){
      adjust_stepsize(dif,(absh/abshlast),neq,k);
      hinvGak = h * invGa[k-1];
      nconhk = 0;
      class_call(new_linearisation(&jac,hinvGak,neq,error_message),
                 error_message,error_message);
      stepstat[4] += 1;
      havrate = _FALSE_;
    }
    /*        Loop for advancing one step */
    nofailed = _TRUE_;
    for( ; ; ){
      gotynew = _FALSE_;    /* is ynew evaluated yet?*/
      while(gotynew==_FALSE_){
        /*Compute the constant terms in the equation for ynew.
          Next FOR lop is just: psi = matmul(dif(:,1:k),(G(1:k) * invGa(k)))*/
        for(ii=1;ii<=neq;ii++){
          psi[ii] = 0.0;
          for(jj=1;jj<=k;jj++){
            psi[ii] += dif[ii][jj]*G[jj-1]*invGa[k-1];
          }
        }
        /* Predict a solution at t+h. */
        tnew = t + h;
        if (done==_TRUE_){
          tnew = tfinal; /*Hit end point exactly. */
        }
        h = tnew - t;          /* Purify h. */
        for(ii=1;ii<=neq;ii++){
          pred[ii] = y[ii];
          for(jj=1;jj<=k;jj++){
            pred[ii] +=dif[ii][jj];
          }
        }
        eqvec(pred,ynew,neq);

        /*The difference, difkp1, between pred and the final accepted
          ynew is equal to the backward difference of ynew of order
          k+1. Initialize to zero for the iteration to compute ynew.
        */

        minnrm = 0.0;
        for(j=1;j<=neq;j++){
          difkp1[j] = 0.0;
          maxtmp = MAX(fabs(ynew[j]),fabs(y[j]));
          invwt[j] = 1.0 / MAX(maxtmp,threshold);
          maxtmp = 100*eps*fabs(ynew[j]*invwt[j]);
          minnrm = MAX(minnrm,maxtmp);
        }
        /* Iterate with simplified Newton method. */
        tooslow = _FALSE_;
        for(iter=1;iter<=maxit;iter++){
          for (ii=1;ii<=neq;ii++){
            tempvec1[ii]=(psi[ii]+difkp1[ii]);
          }
          class_call((*derivs)(tnew,ynew+1,f0+1,parameters_and_workspace_for_derivs,error_message),
                 error_message,error_message);
          stepstat[2] += 1;
          for(j=1;j<=neq;j++){
            rhs[j] = hinvGak*f0[j]-tempvec1[j];
          }

          /*Solve the linear system A*x=del by using the LU decomposition stored in jac.*/
          if (jac.use_sparse){
            funcreturn = sp_lusolve(jac.Numerical, rhs+1, del+1);
            class_test(funcreturn == _FAILURE_,error_message,
            "Failure in sp_lusolve. Possibly singular matrix!");
          }
          else{
            eqvec(rhs,del,neq);
            funcreturn = lubksb(jac.LU,neq,jac.luidx,del);
            class_test(funcreturn == _FAILURE_,error_message,
            "Failure in lubksb. Possibly singular matrix!");
          }

          stepstat[5]+=1;
          newnrm = 0.0;
          for(j=1;j<=neq;j++){
            maxtmp = fabs(del[j]*invwt[j]);
            newnrm = MAX(newnrm,maxtmp);
          }
          for(j=1;j<=neq;j++){
            difkp1[j] += del[j];
            ynew[j] = pred[j] + difkp1[j];
          }
          if (newnrm <= minnrm){
            gotynew = _TRUE_;
            break; /* Break Newton loop */
          }
          else if(iter == 1){
            if (havrate==_TRUE_){
              errit = newnrm * rate / (1.0 - rate);
              if (errit <= 0.05*rtol){
            gotynew = _TRUE_;
            break; /* Break Newton Loop*/
              }
            }
            else {
              rate = 0.0;
            }
          }
          else if(newnrm > 0.9*oldnrm){
            tooslow = _TRUE_;
            break; /*Break Newton lop */
          }
          else{
            rate = MAX(0.9*rate, newnrm / oldnrm);
            havrate = _TRUE_;
            errit = newnrm * rate / (1.0 - rate);
            if (errit <= 0.5*rtol){
              gotynew = _TRUE_;
              break; /* exit newton */
            }
            else if (iter == maxit){
              tooslow = _TRUE_;
              break; /*exit newton */
            }
            else if (0.5*rtol < errit*pow(rate,(maxit-iter))){
              tooslow = _TRUE_;
              break; /*exit Newton */
            }
          }
          oldnrm = newnrm;
        }
        if (tooslow==_TRUE_){
          stepstat[1] += 1;
          /*    ! Speed up the iteration by forming new linearization or reducing h. */
          if (Jcurrent==_FALSE_){
            class_call((*derivs)(t,y+1,f0+1,parameters_and_workspace_for_derivs,error_message),
                       error_message,error_message);
            nfenj=0;
            class_call(numjac((*derivs),t,y,f0,&jac,&nj_ws,abstol,neq,
                       &nfenj,parameters_and_workspace_for_derivs,error_message),
                       error_message,error_message);
            stepstat[3] += 1;
            stepstat[2] += (nfenj + 1);
            Jcurrent = _TRUE_;
          }
          else if (absh <= hmin){
            class_test(absh <= hmin, error_message,
                       "Step size too small: step:%g, minimum:%g, in interval: [%g:%g]\n",
                       absh,hmin,t0,tfinal);
          }
          else{
            abshlast = absh;
            absh = MAX(0.3 * absh, hmin);
            h = tdir * absh;
            done = _FALSE_;
            adjust_stepsize(dif,(absh/abshlast),neq,k);
            hinvGak = h * invGa[k-1];
            nconhk = 0;
          }
          /* A new linearisation is needed in both cases */
          class_call(new_linearisation(&jac,hinvGak,neq,error_message),
                     error_message,error_message);
          stepstat[4] += 1;
          havrate = _FALSE_;
        }
      }
      /*end of while loop for getting ynew
    difkp1 is now the backward difference of ynew of order k+1. */
      err = 0.0;
      for(jj=1;jj<=neq;jj++){
        err = MAX(err,fabs(difkp1[jj]*invwt[jj]));
      }
      err = err * erconst[k-1];
      if (err>rtol){
        /*Step failed */
        stepstat[1]+= 1;
        if (absh <= hmin){
          class_test(absh <= hmin, error_message,
                     "Step size too small: step:%g, minimum:%g, in interval: [%g:%g]\n",
                     absh,hmin,t0,tfinal);
        }
        abshlast = absh;
        if (nofailed==_TRUE_){
          nofailed = _FALSE_;
          hopt = absh * MAX(0.1, 0.833*pow((rtol/err),(1.0/(k+1))));
          if (k > 1){
            errkm1 = 0.0;
            for(jj=1;jj<=neq;jj++){
              errkm1 = MAX(errkm1,fabs((dif[jj][k]+difkp1[jj])*invwt[jj]));
            }
            errkm1 = errkm1*erconst[k-2];
            hkm1 = absh * MAX(0.1, 0.769*pow((rtol/errkm1),(1.0/k)));
            if (hkm1 > hopt){
              hopt = MIN(absh,hkm1);         /* don't allow step size increase */
              k = k - 1;
            }
          }
          absh = MAX(hmin, hopt);
        }
        else{
          absh = MAX(hmin, 0.5 * absh);
        }
        h = tdir * absh;
        if (absh < abshlast){
          done = _FALSE_;
        }
        adjust_stepsize(dif,(absh/abshlast),neq,k);
        hinvGak = h * invGa[k-1];
        nconhk = 0;
        class_call(new_linearisation(&jac,hinvGak,neq,error_message),
                   error_message,error_message);
        stepstat[4] += 1;
        havrate = _FALSE_;
      }
      else {
        break; /* Succesfull step */
      }
    }
    /* End of conditionless FOR loop */
    stepstat[0] += 1;

    /* Update dif: */
    for(jj=1;jj<=neq;jj++){
      dif[jj][k+2] = difkp1[jj] - dif[jj][k+1];
      dif[jj][k+1] = difkp1[jj];
    }
    for(j=k;j>=1;j--){
      for(ii=1;ii<=neq;ii++){
        dif[ii][j] += dif[ii][j+1];
      }
    }
    /** Output **/
    while ((next<tres)&&(tdir * (tnew - t_vec[next]) >= 0.0)){
      /* Do we need to write output? */
      if (tnew==t_vec[next]){
        class_call((*output)(t_vec[next],ynew+1,f0+1,next,parameters_and_workspace_for_derivs,error_message),
                   error_message,error_message);
// MODIFICATION BY LUC
// All print_variables have been moved to the end of time step
/*
    if (print_variables != NULL){
      class_call((*print_variables)(t_vec[next],ynew+1,f0+1,
                    parameters_and_workspace_for_derivs,error_message),
             error_message,error_message);
    }
*/
      }
      else {
        /*Interpolate if we have overshot sample values*/
        interp_from_dif(t_vec[next],tnew,ynew,h,dif,k,yinterp,ypinterp,yppinterp,interpidx,neq,2);

        class_call((*output)(t_vec[next],yinterp+1,ypinterp+1,next,parameters_and_workspace_for_derivs,
                   error_message),error_message,error_message);

      }
      next++;
    }
    /** End of output **/
    if (done==_TRUE_) {
      break;
    }
    klast = k;
    abshlast = absh;
    nconhk = MIN(nconhk+1,maxk+2);
    if (nconhk >= k + 2){
      temp = 1.2*pow((err/rtol),(1.0/(k+1.0)));
      if (temp > 0.1){
        hopt = absh / temp;
      }
      else {
        hopt = 10*absh;
      }
      kopt = k;
      if (k > 1){
        errkm1 = 0.0;
        for(jj=1;jj<=neq;jj++){
          errkm1 = MAX(errkm1,fabs(dif[jj][k]*invwt[jj]));
        }
        errkm1 = errkm1*erconst[k-2];
        temp = 1.3*pow((errkm1/rtol),(1.0/k));
        if (temp > 0.1){
          hkm1 = absh / temp;
        }
        else {
          hkm1 = 10*absh;
        }
        if (hkm1 > hopt){
          hopt = hkm1;
          kopt = k - 1;
        }
      }
      if (k < maxk){
        errkp1 = 0.0;
        for(jj=1;jj<=neq;jj++){
          errkp1 = MAX(errkp1,fabs(dif[jj][k+2]*invwt[jj]));
        }
        errkp1 = errkp1*erconst[k];
        temp = 1.4*pow((errkp1/rtol),(1.0/(k+2.0)));
        if (temp > 0.1){
          hkp1 = absh / temp;
        }
        else {
          hkp1 = 10*absh;
        }
        if (hkp1 > hopt){
          hopt = hkp1;
          kopt = k + 1;
        }
      }
      if (hopt > absh){
        absh = hopt;
        if (k!=kopt){
          k = kopt;
        }
      }
    }
    /* Advance the integration one step. */
    t = tnew;
    eqvec(ynew,y,neq);
    Jcurrent = _FALSE_;

// MODIFICATION BY LUC
    if (print_variables!=NULL){
      class_call((*derivs)(tnew,
                     ynew+1,
                     f0+1,
                     parameters_and_workspace_for_derivs,error_message),
                 error_message,
                 error_message);

        class_call((*print_variables)(tnew,ynew+1,f0+1,
                    parameters_and_workspace_for_derivs,error_message),
                    error_message,error_message);
    }
// end of modification

  }

  /* a last call is compulsory to ensure that all quantitites in
     y,dy,parameters_and_workspace_for_derivs are updated to the
     last point in the covered range */
  class_call((*derivs)(tnew,
                   ynew+1,
                   f0+1,
                   parameters_and_workspace_for_derivs,error_message),
             error_message,
             error_message);

  if (print_variables!=NULL){
    /** If we are printing variables, we must store the final point */
    class_call((*print_variables)(tnew,ynew+1,f0+1,
                  parameters_and_workspace_for_derivs,error_message),
               error_message,error_message);
  }

  if (verbose > 0){
    printf("\n End of evolver. Next=%d, t=%e and tnew=%e.",next,t,tnew);
    printf("\n Statistics: [%d %d %d %d %d %d] \n",stepstat[0],stepstat[1],
       stepstat[2],stepstat[3],stepstat[4],stepstat[5]);
  }

  /** Deallocate memory */

  free(buffer);

  /*     free(f0); */
  /*     free(wt); */
  /*     free(ddfddt); */
  /*     free(pred); */
  /*     free(y); */
  /*     free(invwt); */
  /*     free(rhs); */
  /*     free(psi); */
  /*     free(difkp1); */
  /*     free(del); */
  /*     free(yinterp); */
  /*     free(ypinterp); */
  /*     free(yppinterp); */
  /*     free(tempvec1); */
  /*     free(tempvec2); */

  /*     free(interpidx); */
  /*     free(dif[1]); */
  /*     free(dif); */

  uninitialize_jacobian(&jac);
  uninitialize_numjac_workspace(&nj_ws);
  return _SUCCESS_;

} /*End of program*/

/**********************************************************************/
/* Here are some small routines used in evolver_ndf15:                */
/* "interp_from_dif", "eqvec", "adjust_stepsize", "calc_C",           */
/* "new_linearisation", "relevant_indices", "ludcmp", "lubksb".       */
/**********************************************************************/

void eqvec(double *datavec,double *emptyvec, int n){
  int i;
  for(i=1;i<=n;i++){
    emptyvec[i] = datavec[i];
  }
}

int calc_C(struct jacobian *jac){
  int nz, i, j, k, n, col, row;
  int duplicate;
  int *Ci, *Cp, *Ai, *Ap;
  n = jac->Numerical->n;Ci = jac->Ci;Cp = jac->Cp;
  Ai = jac->spJ->Ai; Ap = jac->spJ->Ap;
  /* Calculate sparse pattern for C = J + J'. We can use jac->Numerical->xi as
     storage, since we can't refactor when jac->repeated_pattern = 0. xi[0..n][0..n].
     At first, Cp[i+1] holds the number of elements in w[i].*/
  for (j=0;j<=n;j++) Cp[j] = 0; /* Clear Cp */
  for (j=0;j<n;j++){
    /* Loop over columns */
    for (i=Ap[j];i<Ap[j+1];i++){
      /* Loop over rows */
      if(Ai[i]!=j){
    /*Don't consider diagonal entries..*/
    /* Add (i,j) to sparsity pattern of C, if it does not exist: */
    duplicate = _FALSE_;
    col = j; row = Ai[i];
    for(k=0;k<Cp[col+1];k++){
      /* Check for duplicates in column col:*/
      if(jac->Numerical->xi[col][k]==row){
        duplicate = _TRUE_;
        break;
      }
    }
    if (!duplicate){
      jac->Numerical->xi[col][Cp[col+1]] = row;
      Cp[col+1]++;
    }
    /* Add (j,i) to sparsity pattern if it does not exist: */
    duplicate = _FALSE_;
    col = Ai[i]; row = j;
    for(k=0;k<Cp[col+1];k++){
      /* Check for dublicates in the Ai[i]'th column: */
      if(jac->Numerical->xi[col][k]==row){
        duplicate = _TRUE_;
        break;
      }
    }
    if (!duplicate){
      jac->Numerical->xi[col][Cp[col+1]] = row;
      Cp[col+1]++;
    }
      }
    }
  }
  /* wi prepared, write sparsity pattern in Ci and Cp:*/
  nz = 0;
  for(j=0;j<n;j++){
    for(k=0;k<Cp[j+1];k++){
      Ci[nz] = jac->Numerical->xi[j][k];
      nz++;
    }
    Cp[j+1] +=Cp[j]; /* We can update Cp here. */
  }
  return _SUCCESS_;
}

/* Subroutine that interpolates from information stored in dif */
int interp_from_difold(double tinterp,double tnew,double *ynew,double h,double **dif,int k, double *yinterp,
            double *ypinterp, double *yppinterp, int* index, int neq, int output){
  /* Output=1: only y_vector. Output=2: y and y prime. Output=3: y, yp and ypp*/
  int i,j,m,l,p,factor;
  double sumj,suml,sump,prodm,s;
  s = (tinterp - tnew)/h;
  if (k==1){
    for(i=1;i<=neq;i++){
      if(index[i]==_TRUE_){
        yinterp[i] = ynew[i] + dif[i][1] * s;
        if (output>1) ypinterp[i] = dif[i][1]/h;
        if (output>2) yppinterp[i] = 0; /* No good estimate can be made of the second derivative */
      }
    }
  }
  else{
    /*This gets tricky */
    for(i=1;i<=neq;i++){
      if(index[i]==_TRUE_){
        /*First the value of the function:    */
        sumj=0.0;
        factor=1;
        for(j=1;j<=k;j++){
          prodm=1.0;
          factor*=j;
          for(m=0;m<j;m++) prodm*=(m+s);
          sumj+=dif[i][j]/factor*prodm;
        }
        yinterp[i] = ynew[i]+sumj;
        /* Now the first derivative: */
        if (output>1){
          factor = 1;
          sumj=0.0;
          for(j=1;j<=k;j++){
            suml = 0.0;
            factor *=j;
            for(l=0;l<j;l++){
              prodm=1.0;
              for(m=0;m<j;m++){
            if(m!=l) prodm*=(m+s);
              }
              suml+=prodm;
            }
            sumj+=dif[i][j]/factor*suml;
          }
          ypinterp[i] = sumj/h;
        }
        /* The second derivative: */
        if (output>2){
          factor=1;
          sumj=0.0;
          for(j=1;j<=k;j++){
            suml=0.0;
            factor*=j;
            for(l=0;l<j;l++){
              sump=0.0;
              for(p=0;p<j;p++){
                if(p!=l){
                  prodm=1.0;
                  for(m=0;m<j;m++){
                    if((m!=l)&&(m!=p)){
                      prodm*=(m+s);
                    }
                  }
                  sump+=prodm;
                }
              }
              suml+=sump;
            }
            sumj+=dif[i][j]/factor*suml;
          }
          yppinterp[i] = sumj/(h*h);
        }
      }
    }
  }
  return _SUCCESS_;
}

/* Subroutine that interpolates from information stored in dif */
int interp_from_dif(double tinterp,
                    double tnew,
                    double *ynew,
                    double h,
                    double **dif,
                    int k,
                    double *yinterp,
                    double *ypinterp,
                    double *yppinterp,
                    int* mask,
                    int neq,
                    int output){
  /* Output=1: only y_vector. Output=2: y and y prime. Output=3: y, yp and ypp*/
  double fact,prod,sumfrac;
  double vecy[5]={0.,0.,0.,0.,0.};
  double vecdy[5]={0.,0.,0.,0.,0.};
  int j, index_x;
  double s, sumtmp, sumtmp2;

  s = (tinterp - tnew)/h;

  prod = 1.0;
  sumfrac = 0.;
  fact = 1.0;
  for (j=0; j<k; j++){
    prod *= (s+j);
    fact *= (j+1);
    sumfrac += 1.0/(s+j);
    vecy[j] = prod/fact;
    vecdy[j] = prod*sumfrac/(h*fact);
  }

  for (index_x=1; index_x<=neq; index_x++){
    if (mask[index_x]==_TRUE_){
      sumtmp = 0;
      sumtmp2 = 0;
      for (j=0; j<k; j++){
        sumtmp += vecy[j]*dif[index_x][j+1];
        sumtmp2 += vecdy[j]*dif[index_x][j+1];
      }
      yinterp[index_x] = ynew[index_x] + sumtmp;
      ypinterp[index_x] = sumtmp2;
    }
  }
  return _SUCCESS_;
}

int adjust_stepsize(double **dif, double abshdivabshlast, int neq,int k){
  double mydifU[5][5]={{-1,-2,-3,-4,-5},{0,1,3,6,10},{0,0,-1,-4,-10},{0,0,0,1,5},{0,0,0,0,-1}};
  double tempvec[5];
  double mydifRU[5][5];
  int ii,jj,kk;

  for(ii=1;ii<=5;ii++) mydifRU[0][ii-1] = -ii*abshdivabshlast;
  for(jj=2;jj<=5;jj++){
    for(ii=1;ii<=5;ii++){
      mydifRU[jj-1][ii-1] = mydifRU[jj-2][ii-1]*(1.0-(1.0+ii*abshdivabshlast)/jj);
    }
  }
  for(ii=0;ii<5;ii++){
    for(kk=0;kk<5;kk++){
      /* Save the i'th row of mydifRU */
      tempvec[kk] = mydifRU[ii][kk];
    }
    for(jj=0;jj<5;jj++){
      /* Now do the matrix multiplication: */
      mydifRU[ii][jj] = 0.0;
      for(kk=0;kk<5;kk++)    mydifRU[ii][jj] += tempvec[kk]*mydifU[kk][jj];
    }
  }

  for(ii=0;ii<neq;ii++){
    for(kk=0;kk<k;kk++){
      /* Save the k first values of the i'th row of dif */
      tempvec[kk] = dif[ii+1][kk+1];
    }
    for(jj=0;jj<k;jj++){
      /* Now do the matrix multiplication: */
      dif[ii+1][jj+1] = 0.0;
      for(kk=0;kk<k;kk++) dif[ii+1][jj+1] += tempvec[kk]*mydifRU[kk][jj];
    }
  }
  return _SUCCESS_;
}

int new_linearisation(struct jacobian *jac,double hinvGak,int neq,ErrorMsg error_message){
  double luparity, *Ax;
  int i,j,*Ap,*Ai,funcreturn;
  if(jac->use_sparse==1){
    Ap = jac->spJ->Ap; Ai = jac->spJ->Ai; Ax = jac->spJ->Ax;
    /* Construct jac->spJ->Ax from jac->xjac, the jacobian:*/
    for(j=0;j<neq;j++){
      for(i=Ap[j];i<Ap[j+1];i++){
        if(Ai[i]==j){
          /* I'm at the diagonal */
          Ax[i] = 1.0-hinvGak*jac->xjac[i];
        }
        else{
          Ax[i] = -hinvGak*jac->xjac[i];
        }
      }
    }
    /* Matrix constructed... */
    if(jac->new_jacobian==_TRUE_){
      /*I have a new pattern, and I have not done a LU decomposition
        since the last jacobian calculation, so    I need to do a full
        sparse LU-decomposition: */
      /* Find the sparsity pattern C = J + J':*/
      calc_C(jac);
      /* Calculate the optimal ordering: */
      sp_amd(jac->Cp, jac->Ci, neq, jac->cnzmax,
         jac->Numerical->q,jac->Numerical->wamd);
      /* if the next line is uncomented, the code uses natural ordering instead of AMD ordering */
      /*jac->Numerical->q = NULL;*/
      funcreturn = sp_ludcmp(jac->Numerical, jac->spJ, 1e-3);
      class_test(funcreturn == _FAILURE_,error_message,
         "Failure in sp_ludcmp. Possibly singular matrix!");
      jac->new_jacobian = _FALSE_;
    }
    else{
      /* I have a repeated pattern, so I can just refactor:*/
      sp_refactor(jac->Numerical, jac->spJ);
    }
  }
  else{
    /* Normal calculation: */
    for(i=1;i<=neq;i++){
      for(j=1;j<=neq;j++){
    jac->LU[i][j] = - hinvGak * jac->dfdy[i][j];
    if(i==j) jac->LU[i][j] +=1.0;
      }
    }
    /*Dense LU decomposition: */
    funcreturn = ludcmp(jac->LU,neq,jac->luidx,&luparity,jac->LUw);
    class_test(funcreturn == _FAILURE_,error_message,
           "Failure in ludcmp. Possibly singular matrix!");
  }
  return _SUCCESS_;
}

/** Helper functions */
int lubksb(double **a, int n, int *indx, double b[]){
  int i,ii=0,ip,j;
  double sum;
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
  return _SUCCESS_;
}


int ludcmp(double **a, int n, int *indx, double *d, double *vv){
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) return _FAILURE_;
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
    big=dum;
    imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
    dum=a[imax][k];
    a[imax][k]=a[j][k];
    a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  return _SUCCESS_;
}

int fzero_Newton(int (*func)(double *x,
                             int x_size,
                             void *param,
                             double *F,
                             ErrorMsg error_message),
                 double *x_inout,
                 double *dxdF,
                 int x_size,
                 double tolx,
                 double tolF,
                 void *param,
                 int *fevals,
                 ErrorMsg error_message){
  /**Given an initial guess x[1..n] for a root in n dimensions,
     take ntrial Newton-Raphson steps to improve the root.
     Stop if the root converges in either summed absolute
     variable increments tolx or summed absolute function values tolf.*/
  int k,i,j,*indx, ntrial=20;
  double errx,errf,d,*F0,*Fdel,**Fjac,*p, *lu_work;
  int has_converged = _FALSE_;
  int funcreturn;
  double toljac = 1e-3;
  double *delx;

  /** All arrays are indexed as [0, n-1] with the exception of p, indx,
      lu_work and Fjac, since they are passed to ludcmp and lubksb. */
  class_alloc(indx, sizeof(int)*(x_size+1), error_message);
  class_alloc(p, sizeof(double)*(x_size+1), error_message);
  class_alloc(lu_work, sizeof(double)*(x_size+1), error_message);
  class_alloc(Fjac, sizeof(double *)*(x_size+1), error_message);
  Fjac[0] = NULL;
  class_alloc(Fjac[1], sizeof(double)*(x_size*x_size+1), error_message);
  for (i=2; i<=x_size; i++){
    Fjac[i] = Fjac[i-1] + x_size;
  }

  class_alloc(F0, sizeof(double)*x_size, error_message);
  class_alloc(delx, sizeof(double)*x_size, error_message);
  class_alloc(Fdel, sizeof(double)*x_size, error_message);

  for (i=1; i<=x_size; i++){
    delx[i-1] = toljac*dxdF[i-1];
  }

  for (k=1;k<=ntrial;k++) {
    /** Compute F(x): */
    /**printf("x = [%f, %f], delx = [%e, %e]\n",
       x_inout[0],x_inout[1],delx[0],delx[1]);*/
    class_call(func(x_inout, x_size, param, F0, error_message),
               error_message, error_message);
    /**    printf("F0 = [%f, %f]\n",F0[0],F0[1]);*/
    *fevals = *fevals + 1;
    errf=0.0; //fvec and Jacobian matrix in fjac.
    for (i=1; i<=x_size; i++)
      errf += fabs(F0[i-1]); //Check function convergence.
    if (errf <= tolF){
      has_converged = _TRUE_;
      break;
    }

    /**
    if (k==1){
      for (i=1; i<=x_size; i++){
        delx[i-1] *= F0[i-1];
      }
    }
    */

    /** Compute the jacobian of F: */
    for (i=1; i<=x_size; i++){
      if (F0[i-1]<0.0)
        delx[i-1] *= -1;
      x_inout[i-1] += delx[i-1];

      /**      printf("x = [%f, %f], delx = [%e, %e]\n",
               x_inout[0],x_inout[1],delx[0],delx[1]);*/
      class_call(func(x_inout, x_size, param, Fdel, error_message),
                 error_message, error_message);
      /**      printf("F = [%f, %f]\n",Fdel[0],Fdel[1]);*/
      for (j=1; j<=x_size; j++)
        Fjac[j][i] = (Fdel[j-1]-F0[j-1])/delx[i-1];
      x_inout[i-1] -= delx[i-1];
    }
    *fevals = *fevals + x_size;

    for (i=1; i<=x_size; i++)
      p[i] = -F0[i-1]; //Right-hand side of linear equations.
    funcreturn = ludcmp(Fjac, x_size, indx, &d, lu_work); //Solve linear equations using LU decomposition.
    class_test(funcreturn == _FAILURE_,error_message,
               "Failure in ludcmp. Possibly singular matrix!");
    funcreturn = lubksb(Fjac, x_size, indx, p);
    class_test(funcreturn == _FAILURE_,error_message,
               "Failure in lubksb. Possibly singular matrix!");
    errx=0.0; //Check root convergence.
    for (i=1; i<=x_size; i++) { //Update solution.
      errx += fabs(p[i]);
      x_inout[i-1] += p[i];
    }
    if (errx <= tolx){
      has_converged = _TRUE_;
      break;
    }
  }

  free(p);
  free(lu_work);
  free(indx);
  free(Fjac[1]);
  free(Fjac);
  free(F0);
  free(delx);
  free(Fdel);

  if (has_converged == _TRUE_){
    return _SUCCESS_;
  }
  else{
    class_stop(error_message, "Newton's method failed to converge. Try improving initial guess on the parameters, decrease the tolerance requirements to Newtons method or increase the precision of the input function.\n");
  }
}

/**********************************************************************/
/* Here are some routines related to the calculation of the jacobian: */
/* "numjac", "initialize_jacobian", "uninitialize_jacobian",                    */
/* "initialize_numjac_workspace", "uninitialize_numjac_workspace".        */
/**********************************************************************/
int numjac(
       int (*derivs)(double x, double * y,double * dy,
             void * parameters_and_workspace, ErrorMsg error_message),
       double t, double *y, double *fval,
       struct jacobian *jac, struct numjac_workspace *nj_ws,
       double thresh, int neq, int *nfe, void * parameters_and_workspace_for_derivs,
       ErrorMsg error_message){
  /*    Routine that computes the jacobian numerically. It is based on the numjac
    implementation in MATLAB, but a feature for recognising sparsity in the
    jacobian and taking advantage of that has been added.
  */
  double eps=1e-16, br=pow(eps,0.875),bl=pow(eps,0.75),bu=pow(eps,0.25);
  double facmin=pow(eps,0.78),facmax=0.1;
  int logjpos, pattern_broken;
  double tmpfac,difmax2=0.,del2,ffscale;
  int i,j,rowmax2;
  double maxval1,maxval2;
  int colmax,group,row,nz,nz2;
  double Fdiff_absrm,Fdiff_new;
  double **dFdy,*fac;
  int *Ap=NULL, *Ai=NULL;

  dFdy = jac->dfdy; /* Assign pointer to dfdy directly for easier notation. */
  fac = jac->jacvec;
  if (jac->use_sparse){
    Ap = jac->spJ->Ap;
    Ai = jac->spJ->Ai;
  }

  /* Set new_jacobian flag: */
  jac->new_jacobian = _TRUE_;

  for(j=1;j<=neq;j++){
    nj_ws->yscale[j] = MAX(fabs(y[j]),thresh);
    nj_ws->del[j] = (y[j] + fac[j] * nj_ws->yscale[j]) - y[j];
  }

  /*Select an increment del for a difference approximation to
    column j of dFdy.  The vector fac accounts for experience
    gained in previous calls to numjac.
  */
  for(j=1;j<=neq;j++){
    if (nj_ws->del[j]==0.0){
      for(;;){
    if (fac[j] < facmax){
      fac[j] = MIN(100*fac[j],facmax);
      nj_ws->del[j] = (y[j] + fac[j]*nj_ws->yscale[j]) - y[j];
      if(nj_ws->del[j]==0.0){
        break;
      }
    }
    else{
      nj_ws->del[j] = thresh;
      break;
    }
      }
    }
  }

  /* keep del pointing into region: */
  for(j=1;j<=neq;j++){
    if (fval[j]>=0.0){
      nj_ws->del[j] = fabs(nj_ws->del[j]);
    }
    else{
      nj_ws->del[j] = -fabs(nj_ws->del[j]);
    }
  }

  /* Sparse calculation?*/
  if ((jac->use_sparse)&&(jac->repeated_pattern >= jac->trust_sparse)){
    /* printf("\n Sparse calculation..neq=%d, has grouping=%d",neq,jac->has_grouping);*/
    /* Everything done sparse'ly. Do we have a grouping? */
    if (jac->has_grouping==0){
      jac->max_group = column_grouping(jac->spJ,jac->col_group,jac->col_wi);
      jac->has_grouping = 1;
    }
    colmax = jac->max_group+1;
    /*    printf("\n                ->groups=%d/%d.",colmax,neq);  */
    for(j=1;j<=colmax;j++){
      /*loop over groups */
      group = j-1;
      for(i=1;i<=neq;i++){
        /*Add y-vector.. */
        nj_ws->ydel_Fdel[i][j] = y[i];
        /*Add del of all groupmembers:*/
        if(jac->col_group[i-1]==group) nj_ws->ydel_Fdel[i][j] +=nj_ws->del[i];
      }
    }
  }
  else{
    /*printf("\n Normal calculation..."); */
    /*Normal calculation: */
    colmax = neq;
    for(j=1;j<=neq;j++){
      for(i=1;i<=neq;i++){
        nj_ws->ydel_Fdel[i][j] = y[i];
      }
      nj_ws->ydel_Fdel[j][j] += nj_ws->del[j];
    }
  }

  /* The next section should work regardless of sparse...*/
  /* Evaluate the function at y+delta vectors:*/
  for(j=1;j<=colmax;j++){
    for(i=1;i<=neq;i++){
      nj_ws->yydel[i] = nj_ws->ydel_Fdel[i][j];
    }
    class_call((*derivs)(t,nj_ws->yydel+1,nj_ws->ffdel+1,
               parameters_and_workspace_for_derivs,error_message),
               error_message,error_message);

    *nfe+=1;
    for(i=1;i<=neq;i++) nj_ws->ydel_Fdel[i][j] = nj_ws->ffdel[i];
  }


  /*Using the Fdel array, form the jacobian and construct max-value arrays.
    First we do it for the sparse case, then for the normal case:*/
  if ((jac->use_sparse)&&(jac->repeated_pattern >= jac->trust_sparse)){
    /* Sparse case:*/
    for(j=0;j<neq;j++){
      /*Loop over columns, and assign corresponding group:*/
      group = jac->col_group[j];
      Fdiff_new = 0.0;
      Fdiff_absrm = 0.0;
      for(i=Ap[j];i<Ap[j+1];i++){
        /* Loop over rows in the sparse matrix */
        row = Ai[i]+1;
        /* Do I want to construct the full jacobian? No, that is ugly..*/
        Fdiff_absrm = MAX(Fdiff_absrm,fabs(Fdiff_new));
        Fdiff_new = nj_ws->ydel_Fdel[row][group+1]-fval[row]; /*Remember to access the column of the corresponding group */
        if (fabs(Fdiff_new)>=Fdiff_absrm){
          nj_ws->Rowmax[j+1] = row;
          nj_ws->Difmax[j+1] = Fdiff_new;
        }
        /* Assign value to sparse rep of jacobian: */
        jac->xjac[i] = Fdiff_new/nj_ws->del[j+1];
      }
      /* The maximum numerical value of Fdel in true column j+1*/
      nj_ws->absFdelRm[j+1] = fabs(nj_ws->ydel_Fdel[nj_ws->Rowmax[j+1]][group+1]);
    }
  }
  else{
    /*Normal case:*/
    for(j=1;j<=neq;j++){
      Fdiff_new = 0.0;
      Fdiff_absrm = 0.0;
      for(i=1;i<=neq;i++){
        Fdiff_absrm = MAX(fabs(Fdiff_new),Fdiff_absrm);
        Fdiff_new = nj_ws->ydel_Fdel[i][j] - fval[i];
        dFdy[i][j] = Fdiff_new/nj_ws->del[j];
        /*Find row maximums:*/
        if(fabs(Fdiff_new)>=Fdiff_absrm){
          /* Found new max location in column */
          nj_ws->Rowmax[j] = i;
          nj_ws->Difmax[j] = fabs(Fdiff_new);
        }
      }
      nj_ws->absFdelRm[j] = fabs(nj_ws->ydel_Fdel[nj_ws->Rowmax[j]][j]);
    }
  }

  /* Adjust fac for next call to numjac. */
  for(i=1;i<=neq;i++){
    nj_ws->absFvalue[i] = fabs(fval[i]);
  }
  for(j=1;j<=neq;j++){
    nj_ws->absFvalueRm[j] = nj_ws->absFvalue[nj_ws->Rowmax[j]];
  }

  logjpos = 0;
  for(j=1;j<=neq;j++){
    if (((nj_ws->absFdelRm[j]<TINY)&&(nj_ws->absFvalueRm[j] < TINY))||(fabs(nj_ws->Difmax[j])<TINY)){
      nj_ws->logj[j] = 1;/*.true.*/
      logjpos = 1;
    }
    else{
      nj_ws->logj[j] = 0;
    }
  }

  if (logjpos ==1){
    for(i=1;i<=neq;i++){
      nj_ws->yydel[i] = y[i];
      nj_ws->Fscale[i] = MAX(nj_ws->absFdelRm[i],nj_ws->absFvalueRm[i]);
    }
    /* If the difference in f values is so small that the column might be just
       ! roundoff error, try a bigger increment. */
    for(j=1;j<=neq;j++){
      if ((nj_ws->logj[j]==1)&&(nj_ws->Difmax[j]<=(br*nj_ws->Fscale[j]))){
        tmpfac = MIN(sqrt(fac[j]),facmax);
        del2 = (y[j] + tmpfac*nj_ws->yscale[j]) - y[j];
        if ((tmpfac!=fac[j])&&(del2!=0.0)){
          if (fval[j] >= 0.0){
            /*! keep del pointing into region */
            del2 = fabs(del2);
          }
          else{
            del2 = -fabs(del2);
          }
          nj_ws->yydel[j] = y[j] + del2;
          class_call((*derivs)(t,nj_ws->yydel+1,nj_ws->ffdel+1,
                     parameters_and_workspace_for_derivs,error_message),
                     error_message,error_message);
          *nfe+=1;
          nj_ws->yydel[j] = y[j];
          rowmax2 = 1;
          Fdiff_new=0.0;
          Fdiff_absrm = 0.0;
          for(i=1;i<=neq;i++){
            Fdiff_absrm = MAX(Fdiff_absrm,fabs(Fdiff_new));
            Fdiff_new = nj_ws->ffdel[i]-fval[i];
            nj_ws->tmp[i] = Fdiff_new/del2;
            if(fabs(Fdiff_new)>=Fdiff_absrm){
              rowmax2 = i;
              difmax2 = fabs(Fdiff_new);
            }
          }
          maxval1 = difmax2*fabs(del2)*tmpfac;
          maxval2 = nj_ws->Difmax[j]*fabs(nj_ws->del[j]);
          if(maxval1>=maxval2){
            /* The new difference is more significant, so
               use the column computed with this increment.
               This depends on wether we are in sparse mode or not: */
            if ((jac->use_sparse)&&(jac->repeated_pattern >= jac->trust_sparse)){
              for(i=Ap[j-1];i<Ap[j];i++) jac->xjac[i]=nj_ws->tmp[Ai[i]+1];
            }
            else{
              for(i=1;i<=neq;i++) dFdy[i][j]=nj_ws->tmp[i];
            }
            /* Adjust fac for the next call to numjac. */
            ffscale = MAX(fabs(nj_ws->ffdel[rowmax2]),nj_ws->absFvalue[rowmax2]);
            if (difmax2 <= bl*ffscale){
              /* The difference is small, so increase the increment. */
              fac[j] = MIN(10*tmpfac, facmax);
            }
            else if(difmax2 > bu*ffscale){
              /* The difference is large, so reduce the increment. */
              fac[j] = MAX(0.1*tmpfac, facmin);
            }
            else{
              fac[j] = tmpfac;
            }
          }
        }
      }
    }
  }
  /* If use_sparse is true but I still don't trust the sparsity pattern, go through the full calculated jacobi-
     matrix, deduce the sparsity pattern, compare with the old pattern, and write the new sparse Jacobi matrix.
     If I do this cleverly, I only have to walk through the jacobian once, and I don't need any local storage.*/

  if ((jac->use_sparse)&&(jac->repeated_pattern < jac->trust_sparse)){
    nz=0; /*Number of non-zeros */
    Ap[0]=0; /*<-Always is.. */
    pattern_broken = _FALSE_;
    for(j=1;j<=neq;j++){
      for(i=1;i<=neq;i++){
        if ((i==j)||(fabs(dFdy[i][j])!=0.0)){
          /* Diagonal or non-zero index found. */
          if (nz>=jac->max_nonzero){
            /* Too many non-zero points to take advantage of sparsity.*/
            jac->use_sparse = 0;
            break;
          }
          /* Test pattern if it is still unbroken: */
          /* Two conditions must be met if the pattern is intact: Ap[j-1]<=nz<Ap[j],
             so that we are in the right column, and (i-1) must exist in column. Ai[nz]*/
          /* We should first test if nz is in the column, otherwise pattern is dead:*/
          if ((pattern_broken==_FALSE_)&&(jac->has_pattern==_TRUE_)){
            if ((nz<Ap[j-1])||(nz>=Ap[j])){
              /* If we are no longer in the right column, pattern is broken for sure. */
              pattern_broken = _TRUE_;
            }
          }
          if ((pattern_broken==_FALSE_)&&(jac->has_pattern==_TRUE_)){
            /* Up to this point, the new jacobian has managed to fit in the old
               sparsity pattern..*/
            if (Ai[nz]!=(i-1)){
              /* The current non-zero rownumber does not fit the current entry in the
             sparse matrix. Pattern MIGHT be broken. Scan ahead in the sparse matrix
             to search for the row entry: (Remember: the indices are sorted..)*/
              pattern_broken = _TRUE_;
              for(nz2=nz; (nz2<Ap[j])&&(Ai[nz2]<=(i-1)); nz2++){
            /* Go through the rest of the column with the added constraint that
               the row index in the sparse matrix should be smaller than the current
               row index i-1:*/
            if (Ai[nz2]==(i-1)){
              /* sparsity pattern recovered.. */
              pattern_broken = _FALSE_;
              nz = nz2;
              break;
            }
            /* Write a zero entry in the sparse matrix, in case we recover pattern. */
            jac->xjac[nz2] = 0.0;
              }
            }
          }
          /* The following works no matter the status of the pattern: */
          /* Write row_number: */
          Ai[nz] = i-1;
          /* Write value: */
          jac->xjac[nz] = dFdy[i][j];
          nz++;
        }
      }
      /* Break this loop too if I have hit max non-zero points: */
      if (jac->use_sparse==_FALSE_) break;
      Ap[j]=nz;
    }
    if (jac->use_sparse==_TRUE_){
      if ((jac->has_pattern==_TRUE_)&&(pattern_broken==_FALSE_)){
        /*New jacobian fitted into the current sparsity pattern:*/
        jac->repeated_pattern++;
        /* printf("\n Found repeated pattern. nz=%d/%d and
           rep.pat=%d.",nz,neq*neq,jac->repeated_pattern); */
      }
      else{
        /*Something has changed (or first run), better still do the full calculation..*/
        jac->repeated_pattern = 0;
      }
      jac->has_pattern = 1;
    }
  }
  return _SUCCESS_;
} /* End of numjac */

int initialize_jacobian(struct jacobian *jac, int neq, ErrorMsg error_message){
  int i;

  if (neq>15){
    jac->use_sparse = 1;
  }
  else{
    jac->use_sparse = 0;
  }
  jac->max_nonzero = (int)(MAX(3*neq,0.20*neq*neq));
  jac->cnzmax = 12*jac->max_nonzero/5;

  /*Maximal number of non-zero entries to be considered sparse */
  jac->repeated_pattern = 0;
  jac->trust_sparse = 4;
  /* Number of times a pattern is repeated before we trust it. */
  jac->has_grouping = 0;
  jac->has_pattern = 0;
  jac->sparse_stuff_initialized=0;

  /*Setup memory for the pointers of the dense method:*/

  class_alloc(jac->dfdy,sizeof(double*)*(neq+1),error_message); /* Allocate vector of pointers to rows of matrix.*/
  class_alloc(jac->dfdy[1],sizeof(double)*(neq*neq+1),error_message);
  jac->dfdy[0] = NULL;
  for(i=2;i<=neq;i++) jac->dfdy[i] = jac->dfdy[i-1]+neq; /* Set row pointers... */

  class_alloc(jac->LU,sizeof(double*)*(neq+1),error_message); /* Allocate vector of pointers to rows of matrix.*/
  class_alloc(jac->LU[1],sizeof(double)*(neq*neq+1),error_message);
  jac->LU[0] = NULL;
  for(i=2;i<=neq;i++) jac->LU[i] = jac->LU[i-1]+neq; /* Set row pointers... */

  class_alloc(jac->LUw,sizeof(double)*(neq+1),error_message);
  class_alloc(jac->jacvec,sizeof(double)*(neq+1),error_message);
  class_alloc(jac->luidx,sizeof(int)*(neq+1),error_message);

  /*Setup memory for the sparse method, if used: */
  if (jac->use_sparse){
    jac->sparse_stuff_initialized = 1;

    jac->xjac=(double*)(jac->luidx+neq+1);
    jac->col_group=(int*)(jac->xjac+jac->max_nonzero);
    jac->col_wi=jac->col_group+neq;
    jac->Cp=jac->col_wi+neq;
    jac->Ci=jac->Cp+neq+1;

    class_alloc(jac->xjac,sizeof(double)*jac->max_nonzero,error_message);
    class_alloc(jac->col_group,sizeof(int)*neq,error_message);
    class_alloc(jac->col_wi,sizeof(int)*neq,error_message);
    class_alloc(jac->Cp,sizeof(int)*(neq+1),error_message);
    class_alloc(jac->Ci,sizeof(int)*jac->cnzmax,error_message);

    class_call(sp_num_alloc(&jac->Numerical, neq,error_message),
           error_message,error_message);

    class_call(sp_mat_alloc(&jac->spJ, neq, neq, jac->max_nonzero,
                error_message),error_message,error_message);

  }

  /* Initialize jacvec to sqrt(eps):*/
  for (i=1;i<=neq;i++) jac->jacvec[i]=1.490116119384765597872e-8;
  return _SUCCESS_;
}

int uninitialize_jacobian(struct jacobian *jac){
  free(jac->dfdy[1]);
  free(jac->dfdy);
  free(jac->LU[1]);
  free(jac->LU);

  free(jac->luidx);
  free(jac->LUw);
  free(jac->jacvec);

  if(jac->sparse_stuff_initialized){
    free(jac->xjac);
    free(jac->col_wi);
    free(jac->col_group);
    free(jac->Cp);
    free(jac->Ci);
    sp_mat_free(jac->spJ);
    sp_num_free(jac->Numerical);
  }
  return _SUCCESS_;
}

int initialize_numjac_workspace(struct numjac_workspace * nj_ws,int neq, ErrorMsg error_message){
  int i,neqp=neq+1;
  /* Allocate vectors and matrices: */

  class_alloc(nj_ws->yscale,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->del,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->Difmax,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->absFdelRm,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->absFvalue,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->absFvalueRm,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->Fscale,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->ffdel,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->yydel,sizeof(double)*neqp,error_message);
  class_alloc(nj_ws->tmp,sizeof(double)*neqp,error_message);

  class_alloc(nj_ws->ydel_Fdel,sizeof(double*)*(neq+1),error_message); /* Allocate vector of pointers to rows of matrix.*/
  class_alloc(nj_ws->ydel_Fdel[1],sizeof(double)*(neq*neq+1),error_message);
  nj_ws->ydel_Fdel[0] = NULL;
  for(i=2;i<=neq;i++) nj_ws->ydel_Fdel[i] = nj_ws->ydel_Fdel[i-1]+neq; /* Set row pointers... */

  class_alloc(nj_ws->logj,sizeof(int)*neqp,error_message);
  class_alloc(nj_ws->Rowmax,sizeof(int)*neqp,error_message);

  /* Done allocating stuff */
  return _SUCCESS_;
}

int uninitialize_numjac_workspace(struct numjac_workspace * nj_ws){
  /* Deallocate vectors and matrices: */
  free(nj_ws->yscale);
  free(nj_ws->del);
  free(nj_ws->Difmax);
  free(nj_ws->absFdelRm);
  free(nj_ws->absFvalue);
  free(nj_ws->absFvalueRm);
  free(nj_ws->Fscale);
  free(nj_ws->ffdel);
  free(nj_ws->yydel);
  free(nj_ws->tmp);

  free(nj_ws->ydel_Fdel[1]);
  free(nj_ws->ydel_Fdel);
  free(nj_ws->logj);
  free(nj_ws->Rowmax);
  return _SUCCESS_;
}
