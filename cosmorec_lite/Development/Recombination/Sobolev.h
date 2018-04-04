//======================================================================================
// Author Jens Chluba July 2007
//======================================================================================

#ifndef SOBOLEV_H
#define SOBOLEV_H

//======================================================================================
double gi_f(const int &l);

//======================================================================================
// exponential factors for recombination routines
//======================================================================================
double B_nu(const double &nu, const double &Tg);
double exp_nu(const double &nu, const double &Tg);
double exp_DE(const double &DE, const double &Tg);

//======================================================================================
// Sobolev optical depth
//======================================================================================
double tau_S(const double &Aij, const double &lij, const int &li, const int &lj, const double &Ni, const double &Nj, const double &H_z);
double dtau_S_dNi(const double &Aij, const double &lij, const int &li, const int &lj, const double &Ni, const double &Nj, const double &H_z);
double dtau_S_dNj(const double &Aij, const double &lij, const int &li, const int &lj, const double &Ni, const double &Nj, const double &H_z);

//======================================================================================
// Sobolev optical depth using statical weights gwi and gwj
//======================================================================================
double tau_S_gw(const double &Aij, const double &lij, const double &gwi, const double &Ni, const double &gwj, const double &Nj, const double &H_z);
double dtau_S_gw_dNi(const double &Aij, const double &lij, const double &gwi, const double &Ni, const double &gwj, const double &Nj, const double &H_z);
double dtau_S_gw_dNj(const double &Aij, const double &lij, const double &gwi, const double &Ni, const double &gwj, const double &Nj, const double &H_z);

//======================================================================================
// Sobolev escape probability
//======================================================================================
double p_ij(const double &tau_S);
double dp_ij_dtau_S(const double &tau_S);

#endif
