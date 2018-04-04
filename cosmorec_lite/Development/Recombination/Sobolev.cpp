//======================================================================================
// Author Jens Chluba July 2007
//======================================================================================

#include <cmath>

#include "Sobolev.h"
#include "routines.h"
#include "physical_consts.h"

using namespace std;

//======================================================================================
double gi_f(const int &l){ return 4.0*(l+0.5); } // == 2(2l+1); factor of 2 due to electron spin

//======================================================================================
// simple routines for the set of equations
//======================================================================================
double B_nu(const double &nu, const double &Tg)
{ return 2.0*const_h*nu*pow(nu/const_cl, 2)/( exp( min(700.0, const_h_kb*nu/Tg) )-1.0); }

double exp_nu(const double &nu, const double &Tg)
{ return exp( min(700.0, const_h_kb*nu/Tg) ); }

double exp_DE(const double &DE, const double &Tg)
{ return exp( min(700.0, DE*1.0e-3/(const_me*Tg)) ); }

//======================================================================================
// Sobolev optical depth
//======================================================================================
double tau_S(const double &Aij, const double &lij, const int &li, const int &lj, const double &Ni, const double &Nj, const double &H_z)
{
    // i is considered to be the 'upper' level (Ni<Nj)
    return Aij*pow(lij, 3)/(8.0*PI*H_z)*fabs(Ni*(Nj/Ni*(li+0.5)/(lj+0.5) - 1.0 ));
}

double dtau_S_dNi(const double &Aij, const double &lij, const int &li, const int &lj, const double &Ni, const double &Nj, const double &H_z)
{
    // i is considered to be the 'upper' level (Ni<Nj)
    return Aij*pow(lij, 3)/(8.0*PI*H_z);
}

double dtau_S_dNj(const double &Aij, const double &lij, const int &li, const int &lj, const double &Ni, const double &Nj, const double &H_z)
{
    // i is considered to be the 'upper' level (Ni<Nj)
    return Aij*pow(lij, 3)/(8.0*PI*H_z)*(li+0.5)/(lj+0.5);
}

//======================================================================================
// Sobolev optical depth using statical weights gwi and gwj
//======================================================================================
double tau_S_gw(const double &Aij, const double &lij, const double &gwi, const double &Ni, const double &gwj, const double &Nj, const double &H_z)
{
    // i is considered to be the 'upper' level (Ni<Nj)
    return Aij*pow(lij, 3)/(8.0*PI*H_z)*fabs(Ni*(Nj/Ni*gwi/gwj - 1.0 ));
}

double dtau_S_gw_dNi(const double &Aij, const double &lij, const double &gwi, const double &Ni, const double &gwj, const double &Nj, const double &H_z)
{
    // i is considered to be the 'upper' level (Ni<Nj)
    return Aij*pow(lij, 3)/(8.0*PI*H_z);
}

double dtau_S_gw_dNj(const double &Aij, const double &lij, const double &gwi, const double &Ni, const double &gwj, const double &Nj, const double &H_z)
{
    // i is considered to be the 'upper' level (Ni<Nj)
    return Aij*pow(lij, 3)/(8.0*PI*H_z)*gwi/gwj;
}

//======================================================================================
// Sobolev escape probability
//======================================================================================
double p_ij(const double &tau_S)
{
    if(tau_S<=1.0e-14) return 1.0-0.5*tau_S;
    if(tau_S>=5.0e+2) return 1.0/tau_S;
    return (1.0-exp(-tau_S))/tau_S;           
}

double dp_ij_dtau_S(const double &tau_S)
{
    if(tau_S<=1.0e-14) return -0.5+tau_S/3.0;
    if(tau_S>=5.0e+2) return -1.0/tau_S/tau_S;
    return -(1.0-exp(-tau_S)*(1.0+tau_S))/tau_S/tau_S;                
}
