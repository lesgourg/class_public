//========================================================================================
// Routines to evaluate the lhs of the coupled system of ODEs for the neutral helium
// recombination. The system is written in the form dXi/dt == g(t, Xi), with Xi=N_i/NHtot
//
// Author: Jens Chluba
// date: July 2010
//========================================================================================
//
#ifndef ODEDEF_XI_HEI_H
#define ODEDEF_XI_HEI_H

#include <vector>

#include "HeI_Atom.h"
#include "Photoionization_cross_section.h"
#include "Cosmos.h"

//==============================================================================
// equivalent of the Ly-a line
//==============================================================================
double Recombination_calc_p_d(double Tg, double A21L, double Ric, Gas_of_HeI_Atoms &HeIA);
double Recombination_calc_p_sc(double Tg, double A21L, double Ric, Gas_of_HeI_Atoms &HeIA);

//==============================================================================
// death-probability for Singlet Ly-n lines
//==============================================================================
double Recombination_calc_p_d_Singlet(int ni, int li, double Tg, double A21L, 
                                      double Ric, Gas_of_HeI_Atoms &HeIA);

double Recombination_calc_p_sc_Singlet(int ni, int li, double Tg, double A21L, 
                                       double Ric, Gas_of_HeI_Atoms &HeIA);

//==============================================================================
// death-probability for any level
//==============================================================================
double Recombination_calc_p_d_level(int ni, int li, int si, int ji, double Tg, 
                                    double A21L, double Ric, Gas_of_HeI_Atoms &HeIA);

double Recombination_calc_p_sc_level(int ni, int li, int si, int ji, double Tg, 
                                     double A21L, double Ric, Gas_of_HeI_Atoms &HeIA);

//==============================================================================
// Jens Chluba, June, 2010
//==============================================================================
double fcorr_Loaded(double Tg);

double DPesc_coh(double eta_S, double eta_c, 
                 double Tg, double Te, 
                 double A21L, double Ric, Gas_of_HeI_Atoms &HeIA, 
                 Photoionization_cross_section_SH &H1s);

double DPesc_coh(double eta_S, double eta_c, 
                 double Tg, double Te, 
                 double pd, Voigtprofile_Dawson &P,
                 Photoionization_cross_section_SH &H1s);

double DPd_coh(double eta_d, double eta_c, 
               double Tg, double Te, 
               Voigtprofile_Dawson &Prof,
               Photoionization_cross_section_SH &H1s);

double DPesc_coh_fcorr(double eta_S, double eta_c, 
                       double Tg, double Te,
                       double A21L, double Ric, Gas_of_HeI_Atoms &HeIA, 
                       Photoionization_cross_section_SH &H1s);

double DPesc_coh_fcorr(double eta_S, double eta_c, 
                       double Tg, double Te, 
                       double pd, Voigtprofile_Dawson &Prof,
                       Photoionization_cross_section_SH &H1s);

//==============================================================================
void evaluate_HI_abs_HeI(const double &z, const double &XH1s, const double &NH, 
                         const double &H_z, const double &Tg, const double &Tm, 
                         Gas_of_HeI_Atoms &HeI_Atoms, Cosmos &cosmos, 
                         Photoionization_cross_section_SH &H1s, 
                         double (*DPesc_func_ptr)(double Tm, double tauS, double eta_c),
                         double &dXe_HeI_HI_abs_dt);

void evaluate_HI_abs_HeI_Intercombination(const double &z, const double &XH1s, 
                                          const double &NH, const double &H_z, 
                                          const double &Tg, const double &Tm, 
                                          Gas_of_HeI_Atoms &HeI_Atoms, 
                                          Cosmos &cosmos, Photoionization_cross_section_SH &H1s, 
                                          double (*DPesc_func_ptr)(double Tm, double tauS, double eta_c),
                                          double &dXe_HeI_HI_abs_dt);

//==============================================================================
// feedback of HeI-lines
//==============================================================================
struct HeLevel_data{
    int n;
    int l;
    int s;
    int j;
    double gw;
    int i;
};

struct feedback_data{
    double lambda21;
    double nu21;
    double A21;
    HeLevel_data upperLevel;
    HeLevel_data lowerLevel;    
};

void evaluate_HeI_feedback(double z, double XH1s, Gas_of_HeI_Atoms &HeI_Atoms, Cosmos &cosmos, vector<feedback_data > &L,
                           double (*DP)(int, double, double, Gas_of_HeI_Atoms&), double *dXi_HeI_dt);

#endif
