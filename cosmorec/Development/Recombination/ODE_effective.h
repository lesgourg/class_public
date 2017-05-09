//======================================================================
// Author Jens Chluba July 2010
// last modification: May 2011
//======================================================================
// May 2011: added simple setup for nD-1s Quadrupole line effect.
//           The approximation of Grin et al. 2010 was used.
//======================================================================
#ifndef ODE_EFFECTIVE_H
#define ODE_EFFECTIVE_H

//======================================================================
// Hydrogen & Helium Atom
//======================================================================
#include "Atom.h"
#include "HeI_Atom.h"

//======================================================================
// beginning of Namespace "ODE_effective"                             //
//======================================================================

namespace ODE_effective 
{
    
    //==================================================================
    // electron temperature equation (rho=Te/Tg)
    //==================================================================
    void evaluate_TM(double z, double Xe, double fHe, double rho, double Tg, double Hz, double &drho_dt);
    //
    void df_drho_evaluate_TM(double z, double Xe, double fHe, double Tg, double Hz, double &drho_dt);
    void df_dXHI1s_evaluate_TM(double z, double Xe, double fHe, double rho, double Tg, double Hz, double &drho_dt);
    void df_dXHeI1s_evaluate_TM(double z, double Xe, double fHe, double rho, double Tg, double Hz, double &drho_dt);
    
    //==================================================================
    // 2s - 1s two-photon decay rate
    //==================================================================
    void evaluate_2s_two_photon_decay(double z, double Tg, double X1s, double X2s, double nu21, 
                                      double A2s1s, double &dX1s_H_dt, double &dX2s_H_dt); 
    //
    void df_dX1s_evaluate_2s_two_photon_decay(double z, double Tg, double nu21, double A2s1s, 
                                              double &dX1s_H_dt, double &dX2s_H_dt); 
    
    void df_dX2s_evaluate_2s_two_photon_decay(double z, double Tg, double nu21, double A2s1s, 
                                              double &dX1s_H_dt, double &dX2s_H_dt); 
    
    //==================================================================
    // nP - 1s channel
    //==================================================================
    void evaluate_Ly_n_channel(double z, double Tg, double X1s, double Xnp, double NH, 
                               double H_z, double A21, double lambda21, 
                               double nu21, double &dX1s_H_dt, 
                               double &dXnp_H_dt, double w=3.0);
    
    void df_dX1s_evaluate_Ly_n_channel(double z, double Tg, double X1s, double Xnp, double NH, 
                                       double H_z, double A21, double lambda21, 
                                       double nu21, double &dX1s_H_dt, 
                                       double &dXnp_H_dt, double w=3.0);
    
    void df_dXnp_evaluate_Ly_n_channel(double z, double Tg, double X1s, double Xnp, double NH, 
                                       double H_z, double A21, double lambda21, 
                                       double nu21, double &dX1s_H_dt, 
                                       double &dXnp_H_dt, double w=3.0);

    //==================================================================
    // nD - 1s quadrupole channel (added May 2011)
    //==================================================================
    void evaluate_nD_1s_Q_channel(double Xnp, double Xnd, double A21, 
                                  double &dXnp_H_dt, 
                                  double &dXnd_H_dt);
    
    void df_dXnp_evaluate_nD_1s_Q_channel(double Xnd, double A21, 
                                          double &dXnp_H_dt, 
                                          double &dXnd_H_dt);
    
    void df_dXnd_evaluate_nD_1s_Q_channel(double Xnp, double A21, 
                                          double &dXnp_H_dt, 
                                          double &dXnd_H_dt);
    
}

//======================================================================
// End of Namespace "ODE_effective"                                   //
//======================================================================




//======================================================================
// beginning of Namespace "ODE_HI_effective"                          //
//======================================================================

namespace ODE_HI_effective 
{
    
    //==================================================================
    // effective recombination and photonionization rates
    //==================================================================
    void evaluate_effective_Rci_Ric_terms(double z, double Xe, double Np, Gas_of_Atoms &HIA, 
                                          vector<double> &Ai, vector<double> &Bi, 
                                          int (*get_HI_index)(int), double *dXi_H_dt);
    //
    void evaluate_effective_Rci_Ric_terms_for_Te_derivative(double z, double Xe, double Nc, 
                                                            vector<double> &Ai, vector<double> &Bi, 
                                                            int (*get_HI_index)(int), double *dXi_H_dt);
    //
    void df_dXHI1s_evaluate_effective_Rci_Ric_terms(double z, double Ne, double Nc, 
                                                    vector<double> &Ai, vector<double> &Bi, 
                                                    int (*get_HI_index)(int), double *dXi_H_dt);
    void df_dXHeI1s_evaluate_effective_Rci_Ric_terms(double z, double Nc, 
                                                     vector<double> &Ai, vector<double> &Bi, 
                                                     int (*get_HI_index)(int), double *dXi_H_dt);
    //
    void df_dXi_evaluate_effective_Rci_Ric_terms(int i, double z, 
                                                 vector<double> &Ai, vector<double> &Bi, 
                                                 int (*get_HI_index)(int), double *dXi_H_dt);
    
    //==================================================================
    // effective Rij-rates
    //==================================================================
    void evaluate_effective_Rij_terms(double z, double Tg, Gas_of_Atoms &HIA, vector<vector<double> > &Rij, int (*get_HI_index)(int), double *dXi_H_dt);
    void df_dXi_evaluate_effective_Rij_terms(int i, double z, double Tg, Gas_of_Atoms &HIA, vector<vector<double> > &Rij, int (*get_HI_index)(int), double *dXi_H_dt);
}

//======================================================================
// End of Namespace "ODE_HI_effective"                                //
//======================================================================




//======================================================================
// beginning of Namespace "ODE_HeI_effective"                         //
//======================================================================

namespace ODE_HeI_effective 
{
    //==================================================================
    // effective recombination and photonionization rates
    //==================================================================
    void evaluate_effective_Rci_Ric_terms(double z, double Xe, double Nc, Gas_of_HeI_Atoms &HeIA, vector<double> &Ai, vector<double> &Bi, int (*get_HeI_index)(int), double *dXi_He_dt);
    //
    void df_dXHI1s_evaluate_effective_Rci_Ric_terms(double z, double Nc, vector<double> &Ai, vector<double> &Bi, int (*get_HeI_index)(int), double *dXi_He_dt);
    void df_dXHeI1s_evaluate_effective_Rci_Ric_terms(double z, double Ne, double Nc, vector<double> &Ai, vector<double> &Bi, int (*get_HeI_index)(int), double *dXi_He_dt);
    //
    void df_dXi_evaluate_effective_Rci_Ric_terms(int i, double z, vector<double> &Ai, vector<double> &Bi, int (*get_HeI_index)(int), double *dXi_He_dt);
    
    //==================================================================
    // effective Rij-rates
    //==================================================================
    void evaluate_effective_Rij_terms(double z, double Tg, Gas_of_HeI_Atoms &HeIA, vector<vector<double> > &Rij, int (*get_HeI_index)(int), double *dXi_He_dt);
    void df_dXi_evaluate_effective_Rij_terms(int i, double z, double Tg, Gas_of_HeI_Atoms &HeIA, vector<vector<double> > &Rij, int (*get_HeI_index)(int), double *dXi_He_dt);
}

//======================================================================
// End of Namespace "ODE_HeI_effective"                               //
//======================================================================

#endif
