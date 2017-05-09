//========================================================================================
// Author: Jens Chluba 
// Date: June 2007
// last modification: Feb 2011
//========================================================================================
#ifndef HEI_ATOM_H
#define HEI_ATOM_H

//========================================================================================
// Here the classes for neutral Helium are defined. It consists of a
// collection of all the levels and their properties
//========================================================================================
#include <vector>
#include "physical_consts.h"
#include "File.h"
#include "Voigtprofiles.h"
#include "Photoionization_cross_section.h"

using namespace std;

//========================================================================================
// to load additional TS lines from n^3P_1 - 1^1S_0 series
//========================================================================================
void Set_flag_Intercombination_lines(int val);

//========================================================================================
// to load n^1D_2-1^1S_0 quadrupole lines for 3<=n<=10 
//========================================================================================
void Set_flag_Quadrupole_lines(int val);

struct Transition_Data_HeI_A                                     
{
  int np, lp, sp, jp;        // transition (nn, ll) --> (np, lp)
  double gwp;                // weight of level 
  //double f;                // oscillator strength
  double Dnu;                // transition frequency
  double DE;                 // transition energy
  double A21;                // Einstein A21 coefficient
  double lambda21;           // wavelength of the transition in cm
  double xxx[5];             // to save some information
};

//========================================================================================
//
// Class Electron_Level_HeI_Singlet
//
//========================================================================================
class Electron_Level_HeI_Singlet
{
private:
// for these S=0
    int nn, ll;                                                // quantum numbers (n, l, s)

    int mess_flag;                                             // for level of message output

    double gw;                                                 // weight of level
    double Dnu;                                                // transition frequency to ground state
    double DE;                                                 // energy difference to ground state
    double nuion;                                              // ionisation frequency in Hz
    double Eion;                                               // ionisation energy in eV
    double Eion_ergs;                                          // ionisation energy in ergs

    //==========================================================
    // storage of additional data for computations 
    //==========================================================
    double Xi;          // Xi=Ni/NH == fraction of atoms in that state
    double Ric;         // photoionization rate

    //==========================================================
    // to store all the dipole transitions and their parameters
    //==========================================================
    Transition_Data_HeI_A ZERO_Data;                    
    vector<Transition_Data_HeI_A> A_values;        // contains all the transition data up to 
                                                   // level nn-1 (i.e. downwards).

    void display_Trans_Data(const Transition_Data_HeI_A &v);

    vector<vector<vector<int> > > Transition_Data_HeI_A_lookup_table;
    void create_Transition_lookup_table();
    
    // weight of level (2S+1)(2L+1)=(2J+1)==(2l+1) for S=0
    double g_l(){return (2.0*ll+1.0); }

public:
    //Konstructors and Destructors
    Electron_Level_HeI_Singlet(){}
    Electron_Level_HeI_Singlet(int n, int l, int mflag=1);   
    ~Electron_Level_HeI_Singlet();
    void init(int n, int l, int mflag=1);        

    int Get_n() const { return nn;} 
    int Get_l() const { return ll;} 
    int Get_S() const { return 0;} 
    int Get_J() const { return ll;} 
    double Get_gw() const {return gw; }

    double Get_DE_1s2() const { return DE;}           // energy difference to ground state in eV
    double Get_Dnu_1s2() const { return Dnu;}         // transition frequency to ground state in Hz
    double Get_E_ion() const { return Eion;}          // ionization energy in eV
    double Get_E_ion_ergs() const {return Eion_ergs;} // ionization energy in ergs
    double Get_nu_ion() const { return nuion;}        // ionization frequency in Hz

    //==============================================
    // storage of additional data for computations 
    //==============================================
    double Get_Xi()  const { return Xi;}             // Xi=Ni/NH == fraction of atoms in that state
    double Get_Ric() const { return Ric;}            // photoionization rate 
    void Set_Xi(double v){ Xi=v; return;}
    void Set_Ric(double v){ Ric=v; return;}

    //==============================================
    // for access to the transition data...
    //==============================================
    void display_downward_transition(int i);
    void display_all_downward_transitions();
    void display_level_information();

    int Get_n_down() const { return A_values.size(); }
    const Transition_Data_HeI_A& Get_Trans_Data(int m) const { return A_values[m]; }
    const Transition_Data_HeI_A& Get_Trans_Data(int np, int lp, int sp, int jp) const;
    void Set_xxx_Transition_Data(int m, int d, double v);

    double Get_A21(int np, int lp, int sp, int jp) const;
    //double Get_f(int np, int lp, int sp, int jp) const;
    double Get_nu21(int np, int lp, int sp, int jp) const;
    double Get_lambda21(int np, int lp, int sp, int jp) const;

    //==============================================
    // Saha-relations with continuum
    //==============================================
    double Ni_NeNc_LTE(double TM)  const;
    double Xi_Saha(double Xe, double Xc, double NH, double TM)  const;
    double Ni_Saha(double Ne, double Nc, double TM)  const; 
};

//========================================================================================
//
// Class Electron_Level_HeI_Triplet
//
//========================================================================================
class Electron_Level_HeI_Triplet
{
private:
// for these S=1
    int nn, ll, jj;                                            // quantum numbers (n, l, s, j)

    int mess_flag;                                             // for level of message output

    double gw;                                                 // weight of level
    double Dnu;                                                // transition frequency to ground state
    double DE;                                                 // energy difference to ground state
    double nuion;                                              // ionisation frequency in Hz
    double Eion;                                               // ionisation energy in eV
    double Eion_ergs;                                          // ionisation energy in ergs

    //==========================================================
    // storage of additional data for computations 
    //==========================================================
    double Xi;          // Xi=Ni/NH == fraction of atoms in that state
    double Ric;         // photoionization rate

    //==========================================================
    // to store all the dipole transitions and their parameters
    //==========================================================
    Transition_Data_HeI_A ZERO_Data;                    
    vector<Transition_Data_HeI_A> A_values;        // contains all the transition data up to 
                                                   // level nn-1 (i.e. downwards).
    
    void display_Trans_Data(const Transition_Data_HeI_A &v);
    
    vector<vector<vector<int> > > Transition_Data_HeI_A_lookup_table;
    void create_Transition_lookup_table();
    
    // weight of level (2S+1)(2L+1)=(2J+1)
    double g_l(){return (2.0*jj+1.0); }

public:
    //Konstructors and Destructors
    Electron_Level_HeI_Triplet(){}
    Electron_Level_HeI_Triplet(int n, int l, int j, int mflag=1);   
    ~Electron_Level_HeI_Triplet();
    void init(int n, int l, int j, int mflag=1);        

    int Get_n() const { return nn;} 
    int Get_l() const { return ll;} 
    int Get_S() const { return 1;} 
    int Get_J() const { return jj;} 
    double Get_gw() const {return gw; }

    double Get_DE_1s2() const { return DE;}           // energy difference to ground state in eV 
    double Get_Dnu_1s2() const { return Dnu;}         // transition frequency to ground state in Hz
    double Get_E_ion() const { return Eion;}         // ionization energy in eV
    double Get_E_ion_ergs() const {return Eion_ergs;}// ionization energy in ergs
    double Get_nu_ion() const { return nuion;}       // ionization frequency in Hz

    //==============================================
    // storage of additional data for computations 
    //==============================================
    double Get_Xi()  const { return Xi;}             // Xi=Ni/NH == fraction of atoms in that state
    double Get_Ric() const { return Ric;}            // photoionization rate 
    void Set_Xi(double v){ Xi=v; return;}
    void Set_Ric(double v){ Ric=v; return;}
                                                   
    //==============================================
    // for access to the transition data...
    //==============================================
    void display_downward_transition(int i);
    void display_all_downward_transitions();
    void display_level_information();

    int Get_n_down() const { return A_values.size();}
    const Transition_Data_HeI_A& Get_Trans_Data(int m) const { return A_values[m]; }
    const Transition_Data_HeI_A& Get_Trans_Data(int np, int lp, int sp, int jp) const;
    void Set_xxx_Transition_Data(int m, int d, double v);
    
    double Get_A21(int np, int lp, int jp, int sp) const;
    //double Get_f(int np, int lp, int jp, int sp) const;
    double Get_nu21(int np, int lp, int sp, int jp) const;
    double Get_lambda21(int np, int lp, int sp, int jp) const;

    //==============================================
    // Saha-relations with continuum
    //==============================================
    double Ni_NeNc_LTE(double TM) const;
    double Xi_Saha(double Xe, double Xc, double NH, double TM) const;
    double Ni_Saha(double Ne, double Nc, double TM) const; 
};

//========================================================================================
//
// Class Electron_Level_HeI_Triplet_no_j
//
//========================================================================================
class Electron_Level_HeI_Triplet_no_j
{
private:
    // not j-resolved
    // for these S=1
    int nn, ll;                                                // quantum numbers (n, l, s)

    int mess_flag;                                             // for level of message output

    double gw;                                                 // weight of level
    double Dnu;                                                // transition frequency to ground state
    double DE;                                                 // energy difference to ground state
    double nuion;                                              // ionisation frequency in Hz
    double Eion;                                               // ionisation energy in eV
    double Eion_ergs;                                          // ionisation energy in ergs

    //==========================================================
    // storage of additional data for computations 
    //==========================================================
    double Xi;          // Xi=Ni/NH == fraction of atoms in that state
    double Ric;         // photoionization rate

    //==========================================================
    // to store all the dipole transitions and their parameters
    //==========================================================
    Transition_Data_HeI_A ZERO_Data;                    
    vector<Transition_Data_HeI_A> A_values;       // contains all the transition data up to 
                                                   // level nn-1 (i.e. downwards).

    void display_Trans_Data(const Transition_Data_HeI_A &v);

    vector<vector<vector<int> > > Transition_Data_HeI_A_lookup_table;
    void create_Transition_lookup_table();

    // weight of level (2S+1)(2L+1)
    double g_l(){return 3.0*(2.0*ll+1.0); }

public:
    //Konstructors and Destructors
    Electron_Level_HeI_Triplet_no_j(){}
    Electron_Level_HeI_Triplet_no_j(int n, int l, int mflag=1);   
    ~Electron_Level_HeI_Triplet_no_j();
    void init(int n, int l, int mflag=1);        

    int Get_n() const { return nn;} 
    int Get_l() const { return ll;} 
    int Get_S() const { return 1;} 
    double Get_gw() const {return gw; }

    double Get_DE_1s2() const { return DE;}           // energy difference to ground state in eV 
    double Get_Dnu_1s2() const { return Dnu;}         // transition frequency to ground state in Hz
    double Get_E_ion() const { return Eion;}         // ionization energy in eV
    double Get_E_ion_ergs() const {return Eion_ergs;}// ionization energy in ergs
    double Get_nu_ion() const { return nuion;}       // ionization frequency in Hz

    //==============================================
    // storage of additional data for computations 
    //==============================================
    double Get_Xi()  const { return Xi;}             // Xi=Ni/NH == fraction of atoms in that state
    double Get_Ric() const { return Ric;}            // photoionization rate 
    void Set_Xi(double v){ Xi=v; return;}
    void Set_Ric(double v){ Ric=v; return;}
                                                   
    //==============================================
    // for access to the transition data...
    //==============================================
    void display_downward_transition(int i);
    void display_all_downward_transitions();
    void display_level_information();

    int Get_n_down() const { return A_values.size();}
    const Transition_Data_HeI_A& Get_Trans_Data(int m) const { return A_values[m]; }
    const Transition_Data_HeI_A& Get_Trans_Data(int np, int lp, int sp, int jp) const;
    void Set_xxx_Transition_Data(int m, int d, double v);
    
    double Get_A21(int np, int lp, int jp, int sp) const;
    //double Get_f(int np, int lp, int jp, int sp) const;
    double Get_nu21(int np, int lp, int sp, int jp) const;
    double Get_lambda21(int np, int lp, int sp, int jp) const;

    //==============================================
    // Saha-relations with continuum
    //==============================================
    double Ni_NeNc_LTE(double TM) const;
    double Xi_Saha(double Xe, double Xc, double NH, double TM) const;
    double Ni_Saha(double Ne, double Nc, double TM) const; 
};

//========================================================================================
//
// Class Atomic_Shell_HeI_Singlet
//
//========================================================================================
class Atomic_Shell_HeI_Singlet
{
private:
    int nn;                                                    // quantum numbers n
    int mess_flag;                                             // for level of message output

    void create_Electron_Levels();

public:
    //Konstructors and Destructors
    Atomic_Shell_HeI_Singlet(){}
    Atomic_Shell_HeI_Singlet(int n, int mflag=1);
    ~Atomic_Shell_HeI_Singlet();
    void init(int n, int mflag=1);
    
    vector<Electron_Level_HeI_Singlet> Angular_Momentum_Level;   // contains electron levels in the shell nn

    int Get_n(){ return nn;} 

    // for access to the transition data...
    void display_all_level();
    void display_level(int i);
    void display_general_data_of_level(int i);
};

//========================================================================================
//
// Class Atomic_Shell_HeI_Triplet
//
//========================================================================================
class Atomic_Shell_HeI_Triplet
{
private:
    int nn;                                                    // quantum numbers n
    int mess_flag;                                             // for level of message output

    void create_Electron_Levels();

public:
    //Konstructors and Destructors
    Atomic_Shell_HeI_Triplet(){}
    Atomic_Shell_HeI_Triplet(int n, int mflag=1);
    ~Atomic_Shell_HeI_Triplet();
    void init(int n, int mflag=1);
    
    vector<vector<Electron_Level_HeI_Triplet> > Angular_Momentum_Level;  // contains electron levels in the shell nn

    int Get_n(){ return nn;} 

    // for access to the transition data...
    void display_all_level();
    void display_level(int i);
    void display_general_data_of_level(int i);
};

//========================================================================================
//
// Class Atomic_Shell_HeI_Triplet_no_j
//
//========================================================================================
class Atomic_Shell_HeI_Triplet_no_j
{
private:
    int nn;                                                    // quantum numbers n
    int mess_flag;                                             // for level of message output

    void create_Electron_Levels();

public:
    //Konstructors and Destructors
    Atomic_Shell_HeI_Triplet_no_j(){}
    Atomic_Shell_HeI_Triplet_no_j(int n, int mflag=1);
    ~Atomic_Shell_HeI_Triplet_no_j();
    void init(int n, int mflag=1);
    
    vector<Electron_Level_HeI_Triplet_no_j > Angular_Momentum_Level;  // contains electron levels in the shell nn

    int Get_n(){ return nn;} 

    // for access to the transition data...
    void display_all_level();
    void display_level(int i);
    void display_general_data_of_level(int i);
};

//========================================================================================
//
// Class Atom_HeI_Singlet
//
//========================================================================================
class Atom_HeI_Singlet
{
private:
 
    int nShells;                                     // number of full shells. Each shell has l=0..n 
    int mess_flag;                                   // for level of message output

    struct Level_nl
    {
    int n;
    int l; 
    };

    vector<Level_nl> Level_Map;     
    void fill_Level_Map();
    void create_Shells();

    //================================================================================
    // to check level access
    //================================================================================
    void check_Level(int n, int l, string mess) const;
    void check_Level(int i, string mess) const;

public:
  
    //Konstructors and Destructors
    Atom_HeI_Singlet();
    Atom_HeI_Singlet(int nS, int mflag=1);
    ~Atom_HeI_Singlet();
    void init(int nS, int mflag=1);
    
    vector<Atomic_Shell_HeI_Singlet> Shell;              // contains the full atomic shells up to nShell
                                                         // comment: Shell[0] is empty!!!

    const Electron_Level_HeI_Singlet& Level(int i) const;
    const Electron_Level_HeI_Singlet& Level(int n, int l) const;
    Electron_Level_HeI_Singlet& Level(int i);
    Electron_Level_HeI_Singlet& Level(int n, int l);
    
    void display_Level_Map();
    int Get_nShells() const { return nShells; }

    int Get_Level_index(int n, int l) const { return l+n*(n-1)/2; } // n>=1
    int Get_n_of_Level(int i) const { return Level_Map[i].n; }
    int Get_l_of_Level(int i) const { return Level_Map[i].l; }
    int Get_s_of_Level(int i) const { return 0; }
    int Get_j_of_Level(int i) const { return Level_Map[i].l; }

    int Get_total_number_of_Levels() const { return Level_Map.size(); }
    int Get_number_of_Levels_until(int nmax) const { return nmax*(nmax+1)/2; }
};

//========================================================================================
//
// Class Atom_HeI_Triplet
//
//========================================================================================
class Atom_HeI_Triplet
{
private:
 
    int nShells;                                     // number of full shells. Each shell has l=0..n 
    int mess_flag;                                            // for level of message output

    struct Level_nl
    {
    int n;
    int l; 
    int j; 
    };

    vector<Level_nl> Level_Map;     
    void fill_Level_Map();
    void create_Shells();

    //================================================================================
    // to check level access
    //================================================================================
    void check_Level(int n, int l, int j, string mess) const;
    void check_Level(int i, string mess) const;

public:
  
    //Konstructors and Destructors
    Atom_HeI_Triplet();
    Atom_HeI_Triplet(int nS, int mflag=1);
    ~Atom_HeI_Triplet();
    void init(int nS, int mflag=1);
    
    vector<Atomic_Shell_HeI_Triplet> Shell;                      // contains the full atomic shells up to nShell
                                                                 // comment: Shell[0] is empty!!!

    const Electron_Level_HeI_Triplet& Level(int i) const;
    const Electron_Level_HeI_Triplet& Level(int n, int l, int j) const;
    Electron_Level_HeI_Triplet& Level(int i);
    Electron_Level_HeI_Triplet& Level(int n, int l, int j);
    
    void display_Level_Map();
    int Get_nShells() const { return nShells; }
    int Get_total_number_of_Levels() const { return Level_Map.size(); }
    int Get_number_of_Levels_until(int nmax) const { return 3*nmax*(nmax+1)/2-2*nmax-1; }
    int Get_Level_index(int n, int l, int j) const; 
    int Get_n_of_Level(int i) const { return Level_Map[i].n; }
    int Get_l_of_Level(int i) const { return Level_Map[i].l; }
    int Get_s_of_Level(int i) const { return 1; }
    int Get_j_of_Level(int i) const { return Level_Map[i].j; }        // n>=2
};

//========================================================================================
//
// Class Atom_HeI_Triplet_no_j
//
//========================================================================================
class Atom_HeI_Triplet_no_j
{
private:
 
    int nShells;                                     // number of full shells. Each shell has l=0..n 
    int mess_flag;                                            // for level of message output

    struct Level_nl
    {
    int n;
    int l; 
    };

    vector<Level_nl> Level_Map;     
    void fill_Level_Map();
    void create_Shells();

    //================================================================================
    // to check level access
    //================================================================================
    void check_Level(int n, int l, string mess) const;
    void check_Level(int i, string mess) const;

public:
  
    //Konstructors and Destructors
    Atom_HeI_Triplet_no_j();
    Atom_HeI_Triplet_no_j(int nS, int mflag=1);
    ~Atom_HeI_Triplet_no_j();
    void init(int nS, int mflag=1);
    
    vector<Atomic_Shell_HeI_Triplet_no_j> Shell;                 // contains the full atomic shells up to nShell
                                                                 // comment: Shell[0] is empty!!!

    const Electron_Level_HeI_Triplet_no_j& Level(int i) const;
    const Electron_Level_HeI_Triplet_no_j& Level(int n, int l) const;
    Electron_Level_HeI_Triplet_no_j& Level(int i);
    Electron_Level_HeI_Triplet_no_j& Level(int n, int l);
    
    void display_Level_Map();
    int Get_nShells() const { return nShells; }
    int Get_total_number_of_Levels() const { return Level_Map.size(); }
    int Get_number_of_Levels_until(int nmax) const;
    int Get_Level_index(int n, int l) const; 
    int Get_n_of_Level(int i) const { return Level_Map[i].n; }
    int Get_l_of_Level(int i) const { return Level_Map[i].l; }
    int Get_s_of_Level(int i) const { return 1; }
};

//========================================================================================
//
// Class Gas_of_HeI_Atoms
//
//========================================================================================
class Gas_of_HeI_Atoms
{
private:
    int indexT;                     // index at which Triplet states start
    int indexT_no_j;                // index at which non-j-resolved Triplet states start
    int nl;                         // total number of levels
    
    void check_transition_data();
    
    //===================================================================================
    Voigtprofile_Dawson phi_HeI_nP_S[11]; 
    Voigtprofile_Dawson phi_HeI_nP_T[11]; 
    Voigtprofile_Dawson phi_HeI_nD_S[11]; 
    //===================================================================================
    
public:
    //===================================================================================
    //Konstructors and Destructors
    //===================================================================================
    Gas_of_HeI_Atoms();
    Gas_of_HeI_Atoms(int nS, int njres, int mflag=1);
    ~Gas_of_HeI_Atoms();
    void init(int nS, int njres, int mflag=1); 
    
    Atom_HeI_Singlet Sing;
    Atom_HeI_Triplet Trip;
    Atom_HeI_Triplet_no_j Trip_no_j;
    
    //===================================================================================
    // to include HeI Lya absorption by hydrogen (added 20.02.2011 by JC)
    //===================================================================================
    Photoionization_cross_section_SH HILyc;
    
    //===================================================================================
    // access to resonance profiles (added 20.02.2011 by JC)
    //===================================================================================
    Voigtprofile_Dawson& nP_S_profile(int n);
    Voigtprofile_Dawson& nP_T_profile(int n);
    Voigtprofile_Dawson& nD_S_profile(int n);
    const Voigtprofile_Dawson& nP_S_profile(int n) const;
    const Voigtprofile_Dawson& nP_T_profile(int n) const;
    const Voigtprofile_Dawson& nD_S_profile(int n) const;

    //===================================================================================
    int Get_total_number_of_Levels() const { return nl;}
    
    int Get_nShells() const { return Sing.Get_nShells(); }
    int Get_indexT() const { return indexT; }
    int Get_indexT_no_j() const { return indexT_no_j; }
    int Get_njres();
    int Get_Level_index(int n, int l, int s, int j) const; 
    
    double Get_A(int n, int l, int s, int j, int np, int lp, int sp, int jp) const;
    //double Get_f(int n, int l, int s, int j, int np, int lp, int sp, int jp) const;
    double Get_nu21(int n, int l, int s, int j, int np, int lp, int sp, int jp) const;
    double Get_lambda21(int n, int l, int s, int j, int np, int lp, int sp, int jp) const;
    
    double Xi(int n, int l, int s, int j) const;
    double Ric(int n, int l, int s, int j) const;
    void Set_Xi(int n, int l, int s, int j, double Xi);
    void Set_Ric(int n, int l, int s, int j, double Ric);
    void Show_NLSJ(int i) const;
    
    int Get_n(int i) const; 
    int Get_l(int i) const;
    int Get_S(int i) const;
    int Get_J(int i) const;
    double Get_gw(int i) const;
    double Get_nu_ion(int i) const;
    
    int Get_n_down() const;
    int Get_n_down(int i) const;
    const Transition_Data_HeI_A& Get_Trans_Data(int i, int m) const;
    const Transition_Data_HeI_A& Get_Trans_Data(int i, int np, int lp, int sp, int jp) const;
    void Set_xxx_Transition_Data(int i, int m, int d, double v);
    
    double Xi(int i) const;
    double Ric(int i) const;
    void Set_Xi(int i, double X);
    void Set_Ric(int i, double X);
    double X_tot() const;
    
    //================================================================================
    // Saha-relations with continuum
    //================================================================================
    double Ni_NeNc_LTE(int i, double TM) const;
    double Xi_Saha(int i, double Xe, double Xc, double NH, double TM) const;
    double Ni_Saha(int i, double Ne, double Nc, double TM) const; 
};

#endif

