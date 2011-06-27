#include "common.h"
#include "quadrature.h"

int main(){
  double I,fxy;
  int i,n;
  double *x,*y,*w;
  ErrorMsg error_message;

  n=2;

  class_call(quadrature_in_rectangle(-2.0,2.0,1.0,4.0,&n,&x,&y,&w,error_message),
	     error_message,
	     error_message);


  for(i=0,I=0.0; i<n; i++){
    // f(x,y) = -x^5+y^8-x^3y^6
    fxy = -pow(x[i],5)+pow(y[i],8)-pow(x[i]*y[i]*y[i],3);
    I += fxy*w[i];
  }

  printf("Integral = %g, exact = %g, relative difference = %g.",
	 I, 116508e0, abs(I-116508)/I);

  free(x);
  free(y);
  free(w);
  return _SUCCESS_;
}


  
