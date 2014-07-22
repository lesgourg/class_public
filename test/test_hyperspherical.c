#include "common.h"
#include "hyperspherical.h"


int main(){

  int sgnK;
  double nu;
  int *lvec;
  int l_size;
  int index_l, index_x, index_nu;
  double xmax, sampling, supersampling;
  double *xvec, *Phi, *nuvec;
  int nu_size;
  double dx;
  double Inum;
  double Iexact, l;
  FILE * outnum, *outexact, *out;
  double maxerr;

  HyperInterpStruct HIS;
  ErrorMsg error_message;

  sgnK = 1;

  switch(sgnK){
  case 0:
    outnum = fopen("I_num.dat","w");
    outexact = fopen("I_exact.dat","w");
  case -1:
    outnum = fopen("I_num_open.dat","w");
    outexact = fopen("I_exact.dat","w");
  case 1:
    outnum = fopen("I_num_closed.dat","w");
    outexact = fopen("I_exact.dat","w");
}


  l_size = 400;
  int l_size_in = l_size;
  nu_size = 100;

  nu = 10.0;

  lvec = malloc(sizeof(int)*l_size);
  out = fopen("lvec.dat","w");
  for (index_l=0; index_l<l_size; index_l++){
    lvec[index_l] = index_l*10;
    fprintf(out,"%d ",lvec[index_l]);
  }
  fclose(out);


  nuvec = malloc(sizeof(double)*nu_size);
  out = fopen("nuvec.dat","w");
  for (index_nu=0; index_nu<nu_size; index_nu++){
    nuvec[index_nu] = (index_nu+3)*35.0;
    fprintf(out,"%.16e ",nuvec[index_nu]);
  }
  fclose(out);

  for (index_nu=0; index_nu<nu_size; index_nu++){

    nu = nuvec[index_nu];

    if (sgnK != 1){
      xmax = 5000./nu;
    }
    else{
      xmax = _PI_*0.499999999;
      int k;
      l_size = l_size_in;
      for (k=0; k<l_size_in; k++){
        if (lvec[k]>=nu){
          l_size = k-1;
          break;
        }
      }
    }


    sampling = 6.0;

    class_call(hyperspherical_HIS_create(sgnK,
                                         nu,
                                         l_size,
                                         lvec,
                                         1e-6,
                                         xmax,
                                         sampling,
                                         lvec[l_size-1]+1,
                                         1e-20,
                                         &HIS,
                                         error_message),
               error_message,
               error_message);

    supersampling = 100*HIS.x_size;
    xvec = malloc(sizeof(double)*supersampling);
    Phi =  malloc(sizeof(double)*supersampling);

    for (index_x=0; index_x<supersampling; index_x++){
      xvec[index_x] = HIS.x[0]+index_x*(HIS.x[HIS.x_size-1]-HIS.x[0])/(supersampling-1.0);
    }


    maxerr = 0.;
    for (index_l=0; index_l<l_size; index_l++){
      class_call(hyperspherical_Hermite6_interpolation_vector_Phi(&HIS,
                                                                  supersampling,
                                                                  index_l,
                                                                  xvec,
                                                                  Phi,
                                                                  error_message),
                 error_message, error_message);

      dx = xvec[1]-xvec[0];
      Inum = 0.0;
      for (index_x=0; index_x<supersampling; index_x++)
        Inum += Phi[index_x];

      Inum *= dx;

      l = lvec[index_l];

      //Iexact = sqrt(_PI_)/(2.*nu)*tgamma(0.5*(l+1.0))/tgamma(1.+0.5*l);
      Iexact = sqrt(_PI_)/(2.*nu)*exp(lgamma(0.5*(l+1.0))-lgamma(1.+0.5*l));

      //    printf("Inum = %.6e, Iexact = %.6e\n",Inum,Iexact);
      fprintf(outnum,"%.16e ",Inum);
      fprintf(outexact,"%.16e ",Iexact);

      maxerr = MAX(maxerr,fabs(Inum-Iexact));

    }

    for (index_l=l_size; index_l<l_size_in; index_l++){
      fprintf(outnum,"%.16e ",0.0);
      fprintf(outexact,"%.16e ",0.0);
    }

    fprintf(outnum,"\n");
    fprintf(outexact,"\n");
    printf("Max relative error: %.16e\n",maxerr);

    hyperspherical_HIS_free(&HIS, error_message);
  }

  fclose(outnum);
  fclose(outexact);

  free(lvec);
  free(xvec);
  free(Phi);

  return _SUCCESS_;
}
