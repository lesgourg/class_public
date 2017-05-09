/*Smits_He_recomb_data.h*/
/*
From 
R.A. Benjamin, E.D. Skillman, and D.P. Smits
ApJ. 514:307-324, 1999
*/


struct Smits_Coeff
{
  int   n;
  int   l;
  int   s;
  double alpha;
  double beta;
};



const struct Smits_Coeff Smits_Coeff_Array_Singlet[15] =
/*  UNITS               cm^3/s      dimensionless   */
/*{index    n,  L,  S,  alpha       beta}       */
{
// singlet
  {1,   0,  0,  1.54E-13,   -0.486},
  {2,   0,  0,  5.55E-15,   -0.451},
  {2,   1,  0,  1.26E-14,   -0.695},
  {3,   0,  0,  1.62E-15,   -0.444},
  {3,   1,  0,  5.01E-15,   -0.700},
  {3,   2,  0,  4.31E-15,   -0.872},
  {4,   0,  0,  7.00E-16,   -0.445},
  {4,   1,  0,  2.43E-15,   -0.708},
  {4,   2,  0,  2.67E-15,   -0.871},
  {4,   3,  0,  1.38E-15,   -1.046},
  {5,   0,  0,  3.66E-16,   -0.454},
  {5,   1,  0,  1.35E-15,   -0.718},
  {5,   2,  0,  1.60E-15,   -0.873},
  {5,   3,  0,  1.22E-15,   -1.048},
  {5,   4,  0,  4.46E-16,   -1.183}    // bug in Porter 4.46E-15 <--> 4.46E-16
};


const struct Smits_Coeff Smits_Coeff_Array_Triplet[14] =
/*  UNITS               cm^3/s      dimensionless   */
/*{index    n,  L,  S,  alpha       beta}       */
{
// triplet
  {2,   0,  1,  1.49E-14,   -0.381},
  {2,   1,  1,  5.61E-14,   -0.639},
  {3,   0,  1,  3.72E-15,   -0.344},
  {3,   1,  1,  1.95E-14,   -0.632},
  {3,   2,  1,  1.33E-14,   -0.868},    // bug in Porter 1.33E-15 <--> 1.33E-14
  {4,   0,  1,  1.50E-15,   -0.328},
  {4,   1,  1,  9.04E-15,   -0.632},
  {4,   2,  1,  8.29E-15,   -0.867},
  {4,   3,  1,  4.15E-15,   -1.046},    // bug in Porter 4.16E-15 <--> 4.15E-15
  {5,   0,  1,  7.51E-16,   -0.327},
  {5,   1,  1,  4.84E-15,   -0.636},
  {5,   2,  1,  4.99E-15,   -0.870},
  {5,   3,  1,  3.65E-15,   -1.048},
  {5,   4,  1,  1.34E-15,   -1.183}
};


double Smits_Rec_Rate(int n, int l, int s, double TK)
{
  int index=n*(n-1)/2+l;
  double T4=TK/1.0e+4;
  if(s==0 && n>=1 && n<=5) return Smits_Coeff_Array_Singlet[index].alpha*pow(T4, Smits_Coeff_Array_Singlet[index].beta);
  if(s==1 && n>=2 && n<=5) return Smits_Coeff_Array_Triplet[index-1].alpha*pow(T4, Smits_Coeff_Array_Triplet[index-1].beta);
  else cout << " no Smits-recombination rate avialable! " << endl;

  return 0.0;
}

double Smits_Rec_Rate_Singlet(int i, double TK)
{
  double T4=TK/1.0e+4;
  int n=Smits_Coeff_Array_Singlet[i].n;

  if(n>=1 && n<=5) return Smits_Coeff_Array_Singlet[i].alpha*pow(T4, Smits_Coeff_Array_Singlet[i].beta);
  else cout << " Smits_Rec_Rate_Singlet: no Smits-recombination rate avialable! " << endl;

  return 0.0;
}

double Smits_Rec_Rate_Triplet(int i, double TK)
{
  double T4=TK/1.0e+4;
  int n=Smits_Coeff_Array_Triplet[i].n;

  if(n>=2 && n<=5) return Smits_Coeff_Array_Triplet[i].alpha*pow(T4, Smits_Coeff_Array_Triplet[i].beta);
  else cout << " Smits_Rec_Rate_Triplet: no Smits-recombination rate avialable! " << endl;

  return 0.0;
}
