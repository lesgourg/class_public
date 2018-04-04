//====================================================================================================================
// Author: Jens Chluba
//
// first implementation: June 2005
// last change         : Aug  2014
//====================================================================================================================
// 01.08.2014: fixed bug for transition data when quadrupole lines are activated
// July 2014: tidied up the code; checked verbosity and recombination rates communication;
// May  2011: Support for electric quadrupole lines were added.
//====================================================================================================================

#ifndef ATOM_H
#define ATOM_H

//====================================================================================================================
// Here the class Atom is defined. It consists of a collection of all the levels and their
// transistion probabilities at a given temperature T. Photoionization cross-sections are 
// included for hydrogenic Atoms
//====================================================================================================================
#include <vector>

#include "Rec_Phot_BB.h"
#include "Voigtprofiles.h"
#include "physical_consts.h"

using namespace std;

//====================================================================================================================
struct Transition_Data_A                                     
{
    unsigned short np, lp, jp; // transition (nn, ll) --> (np, lp)
    double Dnu;                // transition frequency
    double A21;                // Einstein A21 coefficient
    double lambda21;           // wavelength of the transition in cm
};

struct Transition_Data_B                                      
{
    unsigned short int np, lp, jp; // transition (nn, ll) --> (np, lp)
    double DJij;                   // to save spectral distortions
};

//====================================================================================================================
class Electron_Level
{
    private:
        int nn, ll;                                    // quantum numbers (n, l)
        double gw;                                     // statistical weight of level (n, l)
        int nmax;                                      // maximal upper quantum number to calculate transitions
        
        int Z;                                         // charge of the hydrogenic atom
        double Np;                                     // effective number of protons in the nucleus (Np=M_nucleus/mp)
        double mu_red;                                 // reduced mass in units of me
        
        int mess_flag;                                 // for level of message output
        int Recombination_flag;                        // 0: no Recombination rates required
                                                       // 1: Storey & Hummer
                                                       // 2: Karzas & Latter
        
        double Dnu;                                    // transition frequency to ground state
        double DE;                                     // energy difference to ground state
        double nuion;                                  // ionisation frequency in Hz
        double Eion;                                   // ionisation energy in eV
        double Eion_ergs;                              // ionisation energy in ergs
        double A;                                      // total Einstein A coefficient;
        double AE2;                                    // total Einstein A coefficient for quadrupole lines;
        double Gamma;
        double Gamma_Q_E2;
        
        Transition_Data_A ZERO_Data;                   // contains all the transition data up to nn-1 (i.e. downwards)
        vector<Transition_Data_A> A_values_down_lm1;   // The data is ordered with respect to n
        vector<Transition_Data_A> A_values_down_lp1; 
        
        // arrays for quadrupole lines
        vector<Transition_Data_A> A_values_down_lm2;   // The data is ordered with respect to n
        vector<Transition_Data_A> A_values_down_lmp;               
        vector<Transition_Data_A> A_values_down_lp2; 
    
        vector<Transition_Data_B>  B_values;           // contains partial upward transition data up to level nmax.

        void display_Trans_Data(Transition_Data_A &v);
        void display_Trans_Data(Transition_Data_B &v);
        
        void calculate_transition_Data();              // calculates all the dipole transitions up to n=nmax
        void calculate_transition_to_level(int n, int l);
        void calculate_DQ_transition_to_level(int n);
        double calculate_A_tot();
        double calculate_A_tot_E2();
        
        // for conversions between A21, B21, B12, etc... 
        double g_l(int l){return 2.0*(2.0*l+1.0); }
        double f_lu2f_ul(double f_lu, int llow, int lup){ return -g_l(llow)*f_lu/g_l(lup); }
        
        double A21(double f_ul, int nu, int lu, int nl, int ll);
        double B12(double f_lu, int nl, int nu);
        double B21(double f_ul, int nu, int nl);
        
        // for photoionization and recombination rate from/to the level (nn, ll)
        Rec_Phot_BB_SH Interaction_with_Photons_SH;     
        Rec_Phot_BB_SH_QSP Interaction_with_Photons_SH_QSP;     

        bool Atom_activate_Quadrupole_lines;

    public:
        //============================================================================================================
        // Konstructors and Destructors
        //============================================================================================================
        Electron_Level(){}
        Electron_Level(int n, int l, int nm, int Z, double Np, int Rec_flag=0, int mflag=1);   
        ~Electron_Level();
        void init(int n, int l, int nm, int Z, double Np, int Rec_flag=0, int mflag=1);        
        
        //============================================================================================================
        // Konstructor with quadrupole line support
        //============================================================================================================
        Electron_Level(int n, int l, int nm, int Z, double Np, bool Qlines_on, 
                       int Rec_flag=0, int mflag=1);   
    
        void init(int n, int l, int nm, int Z, double Np, bool Qlines_on, 
                  int Rec_flag=0, int mflag=1);        
            
        //============================================================================================================
        double DE_ul(int nu, int nl){return Z*Z*const_EH_inf*mu_red*(1.0/nl/nl-1.0/nu/nu); }           // in eV
        double nu_ul(int nu, int nl){ return Z*Z*const_EH_inf_Hz*mu_red*(1.0/nl/nl-1.0/nu/nu); }       // in Hz
        double E_ion(int n){return Z*Z*const_EH_inf*mu_red/n/n; }                                      // in eV
        double E_ion_ergs(int n){return Z*Z*const_EH_inf_ergs*mu_red/n/n; }                            // in ergs
        double nu_ion(int n){ return Z*Z*const_EH_inf_Hz*mu_red/n/n; }                                 // in Hz
        double lambda_ul(int nu, int nl){ return const_cl/nu_ul(nu, nl); }                             // in cm
        
        int Get_n() const { return nn;} 
        int Get_l() const { return ll;} 
        double Get_gw() const { return gw;} 
        int Get_nmax() const { return nmax;}
        int Get_Z() const{ return Z; } 
        double Get_mu_red() const { return mu_red;}
        double Get_Np() const{ return Np; } 
        double Get_DE_1s() const { return DE;}            // energy difference to ground state 
        double Get_Dnu_1s() const { return Dnu;}          // transition frequency to ground state
        double Get_E_ion() const { return Eion;}          // ionization energy
        double Get_E_ion_ergs() const {return Eion_ergs;} // ionization energy in ergs
        double Get_nu_ion() const { return nuion;}        // ionization frequency
        double Get_A() const { return A;}                 // total down transition rate in 1/sec
        double Get_tau() const { return 1.0/A;}           // life time of level in sec
        double Get_Gamma() const { return Gamma;} 
        double Get_Gamma_E1() const { return Gamma-Gamma_Q_E2;} 
        double Get_Gamma_E2() const { return Gamma_Q_E2;} 
        
        double temp[5];                                   // to save some information
        
        //============================================================================================================
        // functions with fast access to transition data
        //============================================================================================================
        int Get_n_down_all() const { return A_values_down_lm1.size()+A_values_down_lp1.size()
                                           +A_values_down_lm2.size()+A_values_down_lmp.size()
                                           +A_values_down_lp2.size(); }
        
        int Get_n_down_lm1() const { return A_values_down_lm1.size(); }
        int Get_n_down_lp1() const { return A_values_down_lp1.size(); }
        //
        int Get_n_down_lm2() const { return A_values_down_lm2.size(); }
        int Get_n_down_lmp() const { return A_values_down_lmp.size(); }
        int Get_n_down_lp2() const { return A_values_down_lp2.size(); }
        //
        const Transition_Data_A& Get_Data_A_lm1(const int &i) const { return A_values_down_lm1[i]; }
        const Transition_Data_A& Get_Data_A_lp1(const int &i) const { return A_values_down_lp1[i]; }
        //
        const Transition_Data_A& Get_Data_A_lm2(const int &i) const { return A_values_down_lm2[i]; }
        const Transition_Data_A& Get_Data_A_lmp(const int &i) const { return A_values_down_lmp[i]; }
        const Transition_Data_A& Get_Data_A_lp2(const int &i) const { return A_values_down_lp2[i]; }
        //
        const Transition_Data_A& Get_Data_A(const int &np, const int &lp) const;
        void display_all_Trans_Data();
        //
        double Get_A21(const int &np, const int &lp) const;
        double Get_gwp(const int &np, const int &lp) const { return 2.0*(2.0*lp+1.0); }     
        double Get_nu21(const int &np, const int &lp) const;
        double Get_lambda21(const int &np, const int &lp) const;
        double Get_f21(const int &np, const int &lp) const
        { return Get_A21(np, lp)*pow(Get_lambda21(np, lp), 2)/(2.0*FOURPI*const_PIe2_mec); }
        double Get_f12(const int &np, const int &lp) const
        { return Get_A21(np, lp)*pow(Get_lambda21(np, lp), 2)/(2.0*FOURPI*const_PIe2_mec)*(2.0*ll+1.0)/(2.0*lp+1.0); }
        
        int Get_n_up() const { return B_values.size(); }
        int Get_B_data_np(const unsigned int &i){ return B_values[i].np; }
        int Get_B_data_lp(const unsigned int &i){ return B_values[i].lp; }
        int Get_B_data_jp(const unsigned int &i){ return B_values[i].jp; }
        double Get_B_data_DJij(const unsigned int &i){ return B_values[i].DJij; }
        void Set_B_data_DJij(const unsigned int &i, const double &x){ B_values[i].DJij=x; return; }
        void display_all_upward_Trans_Data();

        //============================================================================================================
        // access to photoionization and recombination rates
        //============================================================================================================
        void clear_Interaction_w_photons();
        double Get_nu_ionization();
        
        double R_nl_c(double T_g);                                          // in 1/sec
        double R_c_nl(double T_g, double rho=1.0);                          // in cm^3/sec
        double dR_c_nl_dTe(double T_g, double rho=1.0);                     // in cm^3/sec
        
        //============================================================================================================
        // Gaunt-factor and photoionisation cross section
        //============================================================================================================
        double sig_phot_ion_nuc();
        double sig_phot_ion(double nu);
        double g_phot_ion(double nu);
        double phi_phot_ion(double nu);                                     // ionization profile
        double phi_phot_rec(double nu, double Tg, int ind_flg=1);           // recombination profile
};

//====================================================================================================================
class Atomic_Shell
    {
    private:
        int nn;                                        // quantum numbers n
        int nmax;                                      // maximal upper quantum number to calculate transitions
        int Z;                                         // charge of the hydrogenic atom
        double Np;                                     // effective number of protons in the nucleus (Np=M_nucleus/mp)
        
        int mess_flag;                                 // for level of message output
        
        void create_Electron_Levels(bool Qlines_on, int Rec_flag);
        bool Atom_activate_Quadrupole_lines;
        
    public:
        //============================================================================================================
        // Konstructors and Destructors
        //============================================================================================================
        Atomic_Shell(){}
        Atomic_Shell(int n, int nm, int Z, double Np, int Rec_flag=0, int mflag=1);
        ~Atomic_Shell();
        void init(int n, int nm, int Z, double Np, int Rec_flag=0, int mflag=1);

        //============================================================================================================
        // Konstructor with quadrupole line support
        //============================================================================================================
        Atomic_Shell(int n, int nm, int Z, double Np, bool Qlines_on, int Rec_flag=0, int mflag=1);
        void init(int n, int nm, int Z, double Np, bool Qlines_on, int Rec_flag=0, int mflag=1);
        
        int Get_n(){ return nn;} 
        int Get_nmax(){ return nmax;}
        int Get_Z() const{ return Z; } 
        double Get_Np() const{ return Np; } 
        
                                                       // all the dipole transitions and their properties are stored
        vector<Electron_Level> Angular_Momentum_Level; // contains electron levels in the shell nn
        
        // for access to the transition data...
        void display_all_level();
        void display_level(unsigned int i);
        void display_general_data_of_level(unsigned int i);
    };

//====================================================================================================================
class Atom
    {
    private:
        
        unsigned int nShells;                          // number of full shells. Each shell has l=0..n
        int Z;                                         // charge of the hydrogenic atom
        double Np;                                     // effective number of protons in the nucleus (Np=M_nucleus/mp)
        int mess_flag;                                 // for level of message output
        
        struct Level_nl
        {
            int n;
            int l; 
        };
        
        vector<Level_nl> Level_Map;     
        void fill_Level_Map();
        void create_Shells(bool Qlines_on, int Rec_flag);
        bool Atom_activate_Quadrupole_lines;
       
    public:
        
        //============================================================================================================
        // Konstructors and Destructors
        //============================================================================================================
        Atom();
        Atom(int nS, int Z, double Np, int Rec_flag=0, int mflag=1);
        ~Atom();
        void init(int nS, int Z, double Np, int Rec_flag=0, int mflag=1);

        //============================================================================================================
        // Konstructor with quadrupole line support
        //============================================================================================================
        Atom(int nS, int Z, double Np, bool Qlines_on, int Rec_flag=0, int mflag=1);
        void init(int nS, int Z, double Np, bool Qlines_on, int Rec_flag=0, int mflag=1);
        
        //============================================================================================================
        vector<Atomic_Shell> Shell;                  // contains the full atomic shells up to nShell
                                                     // comment: Shell[0] is empty!!!
        
        const Electron_Level& Level(const unsigned int &i) const;
        const Electron_Level& Level(const unsigned int &n, const unsigned int &l) const;
        Electron_Level& Level(const unsigned int &i);
        Electron_Level& Level(const unsigned int &n, const unsigned int &l);
        
        void display_Level_Map();
        int Get_nShells() const { return nShells; }
        int Get_Z() const{ return Z; } 
        double Get_Np() const{ return Np; } 
        int Get_Level_index(const unsigned int &n, const unsigned int &l) const { return l+n*(n-1)/2; } // n>=1
        int Get_n_of_Level(const unsigned int &i) const { return Level_Map[i].n; }
        int Get_l_of_Level(const unsigned int &i) const { return Level_Map[i].l; }
        unsigned int Get_total_number_of_Levels() const { return Level_Map.size(); }
        int Get_number_of_Levels_until(int nmax) const { return nmax*(nmax+1)/2; }
        
        void display_general_data_of_Shell(unsigned int i);
    };

//====================================================================================================================
class Gas_of_Atoms : public Atom
    {
    private:
        double T;                                             // temperature of the gas
        double N;                                             // number density of nuclei of the specific atom
                                                              // (includes ionized atoms...)
        
        // all the populations of levels X=N_nl/N; 
        // N is the total number of nuclei of the specific atom (includes ionized atoms...)
        vector< vector<double> > Xvec;
        vector< vector<double> > Ric_BB_vec;                  // photoionization rate for level i (T is necessary)
        vector< vector<double> > Rci_BB_vec;                  // recombination rate for level i   (T is necessary)
        vector< vector<double> > dRci_dTe_BB_vec;             // derivative with respect to Te    (T is necessary)
        
        void create_vectors_X(int imax);
        
        Voigtprofile_Dawson phi_HI_Lyn [101];
        Voigtprofile_Dawson phi_HI_nD1s[101];
        bool Atom_activate_Quadrupole_lines;

    public:
        
        Voigtprofile_Dawson& HI_Lyn_profile(int n);
        const Voigtprofile_Dawson& HI_Lyn_profile(int n) const;
        // quadrupole line profiles
        Voigtprofile_Dawson& HI_nD1s_profile(int n);
        const Voigtprofile_Dawson& HI_nD1s_profile(int n) const;
        
        // total cross section for photon interacting with the 1s-state (JC 01.03.2011)
        double compute_sigma_tot_nu_1s(double nu, double Tm);
        // total cross section for photon interacting with shell n (JC 09.03.2011)
        double compute_sigma_tot_nu_shell(int n, double nu, double Tm);

        //============================================================================================================
        // Konstructors and Destructors
        //============================================================================================================
        // Rec_flag: 0 no recombination rates
        // Rec_flag: 1 Storey & Hummer with interpolation (very fast; avoids long recursions for high levels)
        // Rec_flag: 2 Storey & Hummer
        //============================================================================================================
        Gas_of_Atoms();
        Gas_of_Atoms(int nS, int Z, double Np, int Rec_flag=0, int mflag=1);
        ~Gas_of_Atoms();
        void init(int nS, int Z, double Np, int Rec_flag=0, int mflag=1); 
        
        //============================================================================================================
        // Konstructor that adds quadrupole lines
        //============================================================================================================
        Gas_of_Atoms(int nS, int Z, double Np, bool Qlines_on, int Rec_flag=0, int mflag=1);
        void init(int nS, int Z, double Np, bool Qlines_on, int Rec_flag=0, int mflag=1); 
        
        int Get_Level_index(const unsigned int &n, const unsigned int &l) const{ return l+n*(n-1)/2; } // n>=1
        double Get_population_of_level(const unsigned int &i) const;
        void Set_population_of_level(const unsigned int &i, const double &x);
        double Get_population_of_level(const unsigned int &n, const unsigned int &l) const;
        void Set_population_of_level(const unsigned int &n, const unsigned int &l, const double &x);
        
        void Set_T(double TT){ T=TT; return; }
        double Get_T() const { return T; }
        void Set_N(double NN){ N=NN; return; }
        double Get_N() const { return N; }
        
        double Get_A(unsigned int n, unsigned int l, unsigned int np, unsigned int lp) const;
        
        double X(const unsigned int &i){ return Get_population_of_level(i); }
        double X(const unsigned int &n, const unsigned int &l){ return Get_population_of_level(n, l); }
        
        double X_tot() const;
        double X_tot(int nlow) const;                             // sum of all populations with n>=nlow
        double R_rec_tot_BB(int nlow, double T, double rho=1.0);  // total recombination rate from nlow upwards
        double R_phot_tot_BB(int nlow, double T);                 // total photoionization rate from nlow upwards
        
        void clear_Interaction_w_photons();
        // for extensive numerical applications
        void update_Ric_BB(double T);
        void update_Rci_BB(double T, double rho=1.0);
        void update_dRci_dTe_BB(double T, double rho=1.0);
        void update_Ric_BB(int nm, double T);
        void update_Rci_BB(int nm, double T, double rho=1.0);
        void update_dRci_dTe_BB(int nm, double T, double rho=1.0);
        // direct access
        void Set_Ric_rate(const unsigned int &n, const unsigned int &l, const double &Ric); 
        void Set_Rci_rate(const unsigned int &n, const unsigned int &l, const double &Rci); 
        void Set_dRci_dTe_rate(const unsigned int &n, const unsigned int &l, const double &dRci_dTe); 
        
        double Ric_BB(unsigned int i); 
        double Ric_BB(unsigned int n, unsigned int l); 
        double Rci_BB(unsigned int i); 
        double Rci_BB(unsigned int n, unsigned int l); 
        double dRci_dTe_BB(unsigned int i); 
        double dRci_dTe_BB(unsigned int n, unsigned int l); 
        double R_rec_tot_BB(int nlow);            // total recombination rate from nlow upwards   (precalculated)
        double R_phot_tot_BB(int nlow);           // total photoionization rate from nlow upwards (precalculated)
        
        //============================================================================================================
        // Saha-relations with continuum
        //============================================================================================================
        double Ni_NeNc_LTE(unsigned int n, unsigned int l, double TM);
        double Ni_NeNc_LTE(unsigned int i, double TM);
        
        double Xi_Saha(unsigned int n, unsigned int l, double Xe, double Xc, double NH, double TM, double z);
        double Xi_Saha(unsigned int i, double Xe, double Xc, double NH, double TM, double z);
        
        double Ni_Saha(unsigned int n, unsigned int l, double Ne, double Nc, double TM);
        double Ni_Saha(unsigned int i, double Ne, double Nc, double TM);
    };

#endif

//====================================================================================================================
//====================================================================================================================

