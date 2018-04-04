//===================================================================================
// Author: Jens Chluba 
// Date  : July 2007
//===================================================================================
#include <iostream>
#include <string>
#include <cmath>

#include "HeI_Atom.h"
#include "routines.h"
#include "physical_consts.h"
#include "File.h"
#include "Voigtprofiles.h"

using namespace std;

//===================================================================================
const bool check_Transition_Data_Tables=false;

//===================================================================================
// He 1s2 ionization potential in eV
//===================================================================================
const double He1s2_ion=198310.6690*const_cl*const_h/const_e;  

//===================================================================================
// location of the atomic model for helium
//===================================================================================
const string path=COSMORECDIR+"./Development/Helium/Helium.Data.lite/";

//===================================================================================
// Rydberg for hydrogenic levels in neutral helium
//===================================================================================
const double const_EHeI=const_EH_inf/(1.0+const_me_malp);

//===================================================================================
// above this value the triplet states will not be treated j-resolved
// this should not be smaller than 10!!!
//===================================================================================
static int HeI_Atom_njresolved=10;  
static int HeI_Atom_read_transition_data=1;  

//===================================================================================
// files for level energies
//===================================================================================
const string nameE_QD=path+"He4_energies_quantumdefects.dat";

//===================================================================================
// files for transition rates
//===================================================================================
const string nameA_SS=path+"He4_SS.dat";
const string nameA_TS=path+"He4_TS.dat";
const string nameA_TT=path+"He4_TT.dat";
const string nameA_ST=path+"He4_ST.dat";
const string nameA_add=path+"He4_add.dat";

//===================================================================================
// 25th May 2009: added other n^3 P_1 - 1^1 S_0 intecombination line for 3<=n<=10 
//===================================================================================
const string nameA_add_Int=path+"He4_add.intercombination_lines.dat";
int _HeI_add_Intercomb=0;
void Set_flag_Intercombination_lines(int val){ _HeI_add_Intercomb=val; return; }

//===================================================================================
// 30th May 2009: added n^1D_2-1^1S_0 quadrupole lines for 3<=n<=10 
//===================================================================================
const string nameA_add_Quad=path+"He4_Quadrupole_lines.dat";
int _HeI_add_Quad=0;
void Set_flag_Quadrupole_lines(int val){ _HeI_add_Quad=val; return; }


const string name_DM_gap_S=path+"DM_gap_S.dat";
const string name_high_S=path+"high_levels_S.dat";
const string name_DM_gap_T=path+"DM_gap_T.dat";
const string name_high_T=path+"high_levels_T.dat";
const string name_high_T_no_j=path+"high_levels_T.no_j.dat";

//===================================================================================
// Read the tables for transition-rates and energies
//===================================================================================
void skip_header(inputFile &f, int nn)
{
    string str;
    for(int i=0; i<nn; i++) str=f.get_next_line();
    
    return;
}

//===================================================================================
void read_transition_inf(string fname, int n, int l, int s, int j, 
                         vector<Transition_Data_HeI_A> &v)
{
    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    double dum;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            f.get_next();
            dum=f.get_next();
            dd.DE=f.get_next()-dum;
            dd.A21=f.get_next();
            //dd.f=f.get_next();
            
            // set derived values
            dd.Dnu=dd.DE*const_e/const_h;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants
            
            v.push_back(dd);
        }
        
        f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
void read_transition_inf_J(string fname, int n, int l, int s, int j, 
                           vector<Transition_Data_HeI_A> &v, int check=0)
{
    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv, line, dont_write_down_again;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        line=f.get_next_int();

        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            dd.Dnu=f.get_next();
            dd.A21=f.get_next();
            
            dont_write_down_again=0;
            // set derived values
            dd.DE=dd.Dnu*const_h/const_e;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants
            
            if(check==1) 
                for(int m=0; m<(int)v.size(); m++)
                    if(v[m].np==dd.np && v[m].lp==dd.lp && v[m].jp==dd.jp 
                    && v[m].sp==dd.sp && v[m].A21==dd.A21){ dont_write_down_again=1; break; }
            
            if(dont_write_down_again==0) v.push_back(dd);
        }
        
        f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
void read_transition_inf_J_non_j(string fname, int n, int l, int s, int j, 
                                 vector<Transition_Data_HeI_A> &v, int check=0)
{
    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv, line, dont_write_down_again;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        line=f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            dd.Dnu=f.get_next();
            dd.A21=f.get_next();
            
            dont_write_down_again=0;
            // set derived values
            dd.DE=dd.Dnu*const_h/const_e;
            if(dd.jp==-10) dd.gwp=3*(2*dd.lp+1);
            else dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants
            
            if(check==1) 
                for(int m=0; m<(int)v.size(); m++)
                    if(v[m].np==dd.np && v[m].lp==dd.lp && v[m].jp==dd.jp 
                    && v[m].sp==dd.sp && v[m].A21==dd.A21)
                    { dont_write_down_again=1; break; }
            
            if(dont_write_down_again==0) v.push_back(dd);
        }
        
        f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
void read_level_energy_DM(string fname, int n, int l, int s, int j, double &E1)
{
    // open the file with transition information 
    inputFile f(fname);
    
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && lv==l && sv==s && jv==j) 
        {
            // skip next entries
            f.get_next_int();
            f.get_next_int();
            f.get_next_int();
            f.get_next_int();
            
            f.get_next();
            E1=f.get_next();
            f.close();   
            return;
        }
        else 
        {
            // upper level
            nv=f.get_next_int();
            sv=f.get_next_int();
            lv=f.get_next_int();
            jv=f.get_next_int();
            
            if(nv==n && lv==l && sv==s && jv==j) 
            {
                // skip next entries
                f.get_next();
                f.get_next();
                E1=f.get_next();
                f.close();   
                return;
            }
        }
        
        f.get_next_line();
    }
    
    cout << " read_level_energy_DM:\n no data on level: (" 
         << n << ", " << l << ", " << s << ", " << j << ") found in file: " 
         << fname << endl;
    
    E1=-6000.0;
    
    f.close();   
    return;
}

//===================================================================================
void read_level_energy_QD(string fname, int n, int l, int s, int j, double &E1)
{
    // open the file with transition information 
    inputFile f(fname);
    
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    while(!f.iseof())
    {
        // lower level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && lv==l && sv==s && jv==j) 
        {
            E1=f.get_next();  // this energy is in eV and positive relative to the continuum
            E1=He1s2_ion-E1;
            f.close();   
            return;
        }
        
        f.get_next_line();
    }
    
    cout << " read_level_energy_QD:\n no data on level: (" 
         << n << ", " << l << ", " << s << ", " << j << ") found in file: " 
         << fname << endl;
    
    E1=-6000.0;
    
    f.close();   
    return;
}


//===================================================================================
void read_transition_add_Quadrupole(string fname, int n, int l, int s, int j, 
                                    vector<Transition_Data_HeI_A> &v, int mflag)
{
    if(mflag>=1) cout << " read_transition_add_Quadrupole::"
                      << " loading n^1 D_2 - 1^1 S_0 quadrupole transitions. " << endl;
    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    double dum;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            f.get_next();
            dum=f.get_next();
            dd.DE=f.get_next()-dum;
            dd.A21=f.get_next();
            //dd.f=f.get_next();
            double fv=f.get_next();
            
            // set derived values
            dd.Dnu=dd.DE*const_e/const_h;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants

            dd.A21=FOURPI*2.0*const_PIe2_mec*fv/dd.lambda21/dd.lambda21*1.0/5.0;        
            
            cout << " (" << nv << " " << lv << " " << sv << " " << jv << ") --> (" 
                 << dd.np << " " << dd.lp << " " << dd.sp << " " << dd.jp << ") A=" 
                 << dd.A21 << " " << 1.0/(FOURPI*2.0*const_PIe2_mec*1.0e+16) << endl;

            v.push_back(dd);
        }
        
        f.get_next_line();
    }
    
    f.close();   
    return;
}

//===================================================================================
void read_transition_add_TS(string fname, int n, int l, int s, int j, 
                            vector<Transition_Data_HeI_A> &v, int mflag)
{
    if(mflag>=1) cout << " read_transition_add_TS::"
                      << " loading additional n^3 P_1 - 1^1 S_0 transitions. " << endl; 

    inputFile f(fname);
    skip_header(f, 7);
    
    int nv, sv, lv, jv;
    double dum;
    Transition_Data_HeI_A dd;
    for(int k=0; k<5; k++) dd.xxx[k]=0.0;
    
    while(!f.iseof())
    {
        f.get_next_int();
        // lower level
        dd.np=f.get_next_int();
        dd.sp=f.get_next_int();
        dd.lp=f.get_next_int();
        dd.jp=f.get_next_int();
        
        // upper level
        nv=f.get_next_int();
        sv=f.get_next_int();
        lv=f.get_next_int();
        jv=f.get_next_int();
        
        if(nv==n && sv==s && lv==l && jv==j) 
        {
            f.get_next();
            dum=f.get_next();
            dd.DE=f.get_next()-dum;
            dd.A21=f.get_next();
            //dd.f=f.get_next();

            // set derived values
            dd.Dnu=dd.DE*const_e/const_h;
            dd.gwp=(2*dd.jp+1);
            dd.lambda21=const_cl/dd.Dnu; // better use my natural constants
            
            if(mflag>=1) cout << " (" << nv << " " << lv << " " << sv << " " 
                              << jv << ") --> (" << dd.np << " " << dd.lp << " " 
                              << dd.sp << " " << dd.jp 
                              << ") A=" << dd.A21 << endl;

            v.push_back(dd);
        }
        
        f.get_next_line();
    }
    
    f.close();   
    return;
}


//###################################################################################
//
// class: Electron_Level_HeI_Singlet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Electron_Level_HeI_Singlet
//===================================================================================
void Electron_Level_HeI_Singlet::init(int n, int l, int mflag)
{
    nn=n; ll=l; 
    mess_flag=mflag;
    gw=g_l();
    A_values.clear();
    
    //===============================================================================
    // energies
    //===============================================================================
    // quantum defect values
    if(ll<=6 && nn>10 && nn<=30) read_level_energy_QD(nameE_QD, nn, ll, 0, ll, DE);
    // hydrogenic energies for all levels l>=7 & n>=11
    else if((ll>6 && nn>10) || (ll<=6 && nn>30)) DE=He1s2_ion-const_EHeI/nn/nn;
    else if((nn==9 && ll==8) || (nn==10 && ll==7) || (nn==10 && ll==8) || (nn==10 && ll==9)) 
        DE=He1s2_ion-const_EHeI/nn/nn;
    // the rest is in Drake & Morton
    else read_level_energy_DM(nameA_SS, nn, ll, 0, ll, DE);
    //===============================================================================
    
    if(DE==-6000.0){ Dnu=Eion=Eion_ergs=nuion=0.0; cout << " BAAAADDD " << endl; }
    else 
    {
        Dnu=DE*const_e/const_h;
        Eion=He1s2_ion-DE;
        Eion_ergs=Eion*const_e;
        nuion=Eion_ergs/const_h;
        
        //===========================================================================
        // all transitions from that level
        //===========================================================================
        if(HeI_Atom_read_transition_data==1)
        {
            // set Drake and Morton values
            read_transition_inf(nameA_SS, nn, ll, 0, ll, A_values);
            read_transition_inf(nameA_TS, nn, ll, 0, ll, A_values);
            read_transition_inf(nameA_add, nn, ll, 0, ll, A_values);
            //
            if(_HeI_add_Quad>=nn && ll==2) 
                read_transition_add_Quadrupole(nameA_add_Quad, nn, ll, 0, ll, A_values, mflag); 
            // hydrogenic
            read_transition_inf_J(name_DM_gap_S, nn, ll, 0, ll, A_values);
            if(nn>10) read_transition_inf_J(name_high_S, nn, ll, 0, ll, A_values);
        }
    }
    
    create_Transition_lookup_table();
    
    if(mess_flag>0)
    {
        cout << "\n %##############################################################%" << endl;
        
        display_level_information();
        cout << endl;
        display_all_downward_transitions();
        
        cout << " %##############################################################%" << endl;
    }
    
    return;
}

//===================================================================================
Electron_Level_HeI_Singlet::Electron_Level_HeI_Singlet(int n, int l, int mflag)
{ init(n, l, mflag); }

//===================================================================================
Electron_Level_HeI_Singlet::~Electron_Level_HeI_Singlet()
{ A_values.clear(); }

//===================================================================================
void Electron_Level_HeI_Singlet::create_Transition_lookup_table()
{
    ZERO_Data.np=ZERO_Data.sp=ZERO_Data.lp=ZERO_Data.jp=ZERO_Data.gwp=0;
    ZERO_Data.A21=ZERO_Data.lambda21=ZERO_Data.Dnu=0.0;
    
    Transition_Data_HeI_A_lookup_table.clear();
    vector<int> dum1;
    vector<int> dum2;
    vector<vector<int> > dum22;

    //=============================================================
    // get all transitions to s=0
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==0) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    int nmin=100000, nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;

    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);

    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
        
    //=============================================================
    // get all transitions to s=1
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==1) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    nmin=100000; nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    // check the data
    if(check_Transition_Data_Tables)
    {
        Transition_Data_HeI_A T;
        for(int m=0; m<Get_n_down(); m++)
        {
            T=Get_Trans_Data(A_values[m].np, A_values[m].lp, A_values[m].sp, A_values[m].jp);
//          cout << m << " " << T.A21 << " " << A_values[m].A21 << endl;
            if(T.np==0 || T.A21-A_values[m].A21!=0.0) wait_f_r("create_Transition_lookup_table::NOOOO");        
            if(T.np>nn){ cout << "create_Transition_lookup_table::Singlet:"
                              << " there is some transtions with np>nn. Here np= " 
                              << T.np << " and nn= " << nn << endl; wait_f_r(); }
        }
    }
    
    return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::display_Trans_Data(const Transition_Data_HeI_A &v)
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Singlet::display_Trans_Data:\n % (" 
       << nn << ", " << ll << ", 0, " << ll << ") --> (" 
       << v.np << ", " << v.lp << ", " << v.sp << ", " << v.jp << "), gw:" << gw << endl;
  cout << " %==============================================================%" << endl;
  cout << scientific << " Dnu: " << v.Dnu*1.0e-9 << " GHz, DE: " << v.DE << " eV, lambda21: " 
       << v.lambda21 << " cm"<< endl;
  //cout << " A21: " << v.A21  << " sec^-1, f: " << v.f << endl << endl; 
  cout << scientific << " A21: " << v.A21  << " sec^-1, gw: " << v.gwp << endl << endl; 
  
  return;
}

void Electron_Level_HeI_Singlet::display_level_information()
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Singlet::display_level_information:\n % (n, l, s, j) == (" 
       << nn << ", " << ll << ", 0, " << ll << "), gw:" << gw  
       << endl;
  cout << " %==============================================================%" << endl;
  
  cout << "\n gw: " << gw << endl;
  cout << scientific << " Dnu_1s: " << Dnu*1.0e-9 << " GHz, nuion: " << nuion*1.0e-9 << " GHz" << endl;
  cout << scientific << " DE_1s: " << DE << " eV, Eion: " << Eion << " eV " << endl;
  cout << " Number of transitions: " << Get_n_down() << endl << endl;
  cout << " %==============================================================%" << endl;

  return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::display_downward_transition(int i)
{
    if(i>=(int)A_values.size())
    { 
    cout << " Electron_Level::display_downward_transition::no transition tabulated " << i << endl;
    return;
    }
    
    display_Trans_Data(A_values[i]);

    return;
}

//===================================================================================
void Electron_Level_HeI_Singlet::display_all_downward_transitions()
{
    for(int l=0; l<(int)A_values.size(); l++) display_Trans_Data(A_values[l]);

    return;
}

//======================================================================
const Transition_Data_HeI_A& Electron_Level_HeI_Singlet::Get_Trans_Data(int np, int lp, 
                                                                        int sp, int jp) const
{
    if(sp==0 && (int)Transition_Data_HeI_A_lookup_table[0].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[0][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[0][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[0][np][m]];      
            }
        }
    }
    if(sp==1 && (int)Transition_Data_HeI_A_lookup_table[1].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[1][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[1][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[1][np][m]];      
            }
        }
    }
    
    return ZERO_Data;
}

//===================================================================================
double Electron_Level_HeI_Singlet::Get_A21(int np, int lp, int sp, int jp) const
{
  for(int i=0; i<(int)A_values.size(); i++) 
    {
      if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
          return A_values[i].A21;
      //else return 0.0;
    }

  return 0.0;
}

/*
double Electron_Level_HeI_Singlet::Get_f(int np, int lp, int sp, int jp) const
{
  for(int i=0; i<(int)A_values.size(); i++) 
    {
      if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
          return A_values[i].f;
      else return 0.0;
    }

  return 0.0;
}
*/

//===================================================================================
double Electron_Level_HeI_Singlet::Get_nu21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].Dnu;
        //else return 0.0;
    }
    
    return 0.0;
}

//===================================================================================
double Electron_Level_HeI_Singlet::Get_lambda21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
//      cout << i << " " 
//           << A_values[i].np << " " << A_values[i].lp << " " 
//           << A_values[i].sp << " " << A_values[i].jp << endl;
        
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].lambda21;
        //else return 0.0;
    }
    
    return 0.0;
}

//===================================================================================
void Electron_Level_HeI_Singlet::Set_xxx_Transition_Data(int m, int d, double v)
{ 
  if(d<=5)
    {
      A_values[m].xxx[d]=v; 
    }
  else cout << " Electron_Level_HeI_Singlet::Set_xxx_Transition_Data: check index " << endl;

  return; 
}


//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
double Electron_Level_HeI_Singlet::Ni_NeNc_LTE(double TM) const
{ 
  // HeI --> gc=2 --> ge*gc=4  
  return gw/4.0*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*TM, -1.5)
               *exp(Eion_ergs/const_kB/TM );
}

double Electron_Level_HeI_Singlet::Xi_Saha(double Xe, double Xc, double NH, double TM) const
{ return Xe*Xc*NH*Ni_NeNc_LTE(TM); }

double Electron_Level_HeI_Singlet::Ni_Saha(double Ne, double Nc, double TM) const
{ return Ne*Nc*Ni_NeNc_LTE(TM); }


//###################################################################################
//
// class: Electron_Level_HeI_Triplet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Electron_Level_HeI_Triplet
//===================================================================================
void Electron_Level_HeI_Triplet::init(int n, int l, int j, int mflag)
{
    nn=n; ll=l; jj=j; 
    mess_flag=mflag;
    gw=g_l();
    A_values.clear();
    
    //===============================================================================
    // energies
    //===============================================================================
    // quantum defect values
    if(ll<=6 && nn>10 && nn<=30) read_level_energy_QD(nameE_QD, nn, ll, 1, jj, DE);
    // hydrogenic energies
    else if((ll>6 && nn>10) || (ll<=6 && nn>30)) DE=He1s2_ion-const_EHeI/nn/nn;
    else if((nn==8 && ll==7 && jj==6)) read_level_energy_DM(nameA_TS, nn, ll, 1, jj, DE);
    else if((nn==8 && ll==7) || (nn==9 && ll==8) || 
            (nn==10 && ll==7)|| (nn==10 && ll==8)|| (nn==10 && ll==9)) DE=He1s2_ion-const_EHeI/nn/nn;
    //
    // the rest is in Drake & Morton
    else read_level_energy_DM(nameA_TT, nn, ll, 1, jj, DE);
    //===============================================================================
    
    if(DE==-6000.0){ Dnu=Eion=Eion_ergs=nuion=0.0; cout << " BAAAADDD " << endl; }
    else
    {
        Dnu=DE*const_e/const_h;
        Eion=He1s2_ion-DE;
        Eion_ergs=Eion*const_e;
        nuion=Eion_ergs/const_h;
        
        //===========================================================================
        // all transitions from that level
        //===========================================================================
        if(HeI_Atom_read_transition_data==1)
        {
            // set Drake and Morton values (j-resolved)
            read_transition_inf(nameA_TT, nn, ll, 1, jj, A_values);
            read_transition_inf(nameA_ST, nn, ll, 1, jj, A_values);
            read_transition_inf(nameA_add, nn, ll, 1, jj, A_values);
            if(_HeI_add_Intercomb>=nn && nn>2 && ll==1 && jj==1) 
                read_transition_add_TS(nameA_add_Int, nn, ll, 1, jj, A_values, mflag); 
            // hydrogenic for the gap in Drake and Morton (j-resolved)
            read_transition_inf_J(name_DM_gap_T, nn, ll, 1, jj, A_values);
            // hydrogenic
            if(nn>10) read_transition_inf_J(name_high_T, nn, ll, 1, jj, A_values);
        }
    }
    
    create_Transition_lookup_table();

    if(mess_flag>0)
    {
        cout << "\n %##############################################################%" << endl;
        
        display_level_information();
        cout << endl;
        display_all_downward_transitions();
        
        cout << " %##############################################################%" << endl;
    }
    
    return;
}

//===================================================================================
Electron_Level_HeI_Triplet::Electron_Level_HeI_Triplet(int n, int l, int j, int mflag)
{ init(n, l, j, mflag); }

Electron_Level_HeI_Triplet::~Electron_Level_HeI_Triplet()
{ A_values.clear(); }

void Electron_Level_HeI_Triplet::create_Transition_lookup_table()
{
    ZERO_Data.np=ZERO_Data.sp=ZERO_Data.lp=ZERO_Data.jp=ZERO_Data.gwp=0;
    ZERO_Data.A21=ZERO_Data.lambda21=ZERO_Data.Dnu=0.0;
    
    Transition_Data_HeI_A_lookup_table.clear();
    vector<int> dum1;
    vector<int> dum2;
    vector<vector<int> > dum22;
    
    //=============================================================
    // get all transitions to s=0
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==0) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    int nmin=100000, nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    //=============================================================
    // get all transitions to s=1
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==1) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    nmin=100000; nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    // check the data
    if(check_Transition_Data_Tables)
    {
        Transition_Data_HeI_A T;
        for(int m=0; m<Get_n_down(); m++)
        {
            T=Get_Trans_Data(A_values[m].np, A_values[m].lp, A_values[m].sp, A_values[m].jp);
//          cout << m << " " << T.A21 << " " << A_values[m].A21 << " " << T.np << " " 
//               << T.lp << " " << T.sp << " " << T.jp << " <-- " 
//               << nn << " " << ll << " " << 1 << " " << jj << endl;
            
            if(T.np==0 || T.A21-A_values[m].A21!=0.0) 
                wait_f_r(" create_Transition_lookup_table::NOOOO");
            if(T.np>nn)
            { 
                cout << " create_Transition_lookup_table::Triplet:"
                     << " there is some transtions with np>nn. Here np= " << T.np 
                     << " and nn= " << nn << endl; 
                wait_f_r(); 
            }
        }
    }
    
    return;
}

//===================================================================================
void Electron_Level_HeI_Triplet::display_Trans_Data(const Transition_Data_HeI_A &v)
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet::display_Trans_Data:\n % (" 
       << nn << ", " << ll << ", 1, " << jj << ") --> (" 
       << v.np << ", " << v.lp << ", " << v.sp << ", " << v.jp << "), gw:" << gw << endl;
  cout << " %==============================================================%" << endl;
  cout << scientific << " Dnu: " << v.Dnu*1.0e-9 << " GHz, DE: " << v.DE << " eV, lambda21: " 
       << v.lambda21 << " cm"<< endl;
  //cout << " A21: " << v.A21  << " sec^-1, f: " << v.f << endl << endl; 
  cout << scientific << " A21: " << v.A21  << " sec^-1, gw: " << v.gwp << endl << endl; 
  
  return;
}

//===================================================================================
void Electron_Level_HeI_Triplet::display_level_information()
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet::display_level_information:\n % (n, l, s, j) == (" 
       << nn << ", " << ll << ", 1, " << jj << "), gw:" << gw  
       << endl;
  cout << " %==============================================================%" << endl;
  
  cout << "\n gw: " << gw << endl;
  cout << scientific << " Dnu_1s: " << Dnu*1.0e-9 << " GHz, nuion: " 
       << nuion*1.0e-9 << " GHz" << endl;
  cout << scientific << " DE_1s: " << DE << " eV, Eion: " << Eion << " eV " << endl;
  cout << " Number of transitions: " << Get_n_down() << endl << endl;
  cout << " %==============================================================%" << endl;

  return;
}

void Electron_Level_HeI_Triplet::display_downward_transition(int i)
{
    if(i>=(int)A_values.size())
    { 
        cout << " Electron_Level_HeI_Triplet::display_downward_transition::"
             << " no transition tabulated " << i << endl;
        return;
    }
    
    display_Trans_Data(A_values[i]);

    return;
}

void Electron_Level_HeI_Triplet::display_all_downward_transitions()
{
    for(int l=0; l<(int)A_values.size(); l++) display_Trans_Data(A_values[l]);

    return;
}

const Transition_Data_HeI_A& Electron_Level_HeI_Triplet::Get_Trans_Data(int np, int lp, 
                                                                        int sp, int jp) const
{
    if(sp==0 && (int)Transition_Data_HeI_A_lookup_table[0].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[0][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {           
            if(Transition_Data_HeI_A_lookup_table[0][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[0][np][m]];      
            }
        }
    }
    if(sp==1 && (int)Transition_Data_HeI_A_lookup_table[1].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[1][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[1][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[1][np][m]];      
            }
        }
    }
    
    return ZERO_Data;
}

//======================================================================
double Electron_Level_HeI_Triplet::Get_A21(int np, int lp, int sp, int jp) const
{
  for(int i=0; i<(int)A_values.size(); i++) 
    {
      if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
          return A_values[i].A21;
      //else return 0.0;
    }

  return 0.0;
}

/*
double Electron_Level_HeI_Triplet::Get_f(int np, int lp, int sp, int jp) const
{
  for(int i=0; i<(int)A_values.size(); i++) 
    {
      if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
          return A_values[i].f;
      else return 0.0;
    }

  return 0.0;
}
*/

double Electron_Level_HeI_Triplet::Get_nu21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].Dnu;
        //else return 0.0;
    }
    
    return 0.0;
}

double Electron_Level_HeI_Triplet::Get_lambda21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].lambda21;
        //else return 0.0;
    }
    
    return 0.0;
}

void Electron_Level_HeI_Triplet::Set_xxx_Transition_Data(int m, int d, double v) 
{ 
  if(d<=5)
    {
      A_values[m].xxx[d]=v; 
    }
  else cout << " Electron_Level_HeI_Triplet::Set_xxx_Transition_Data: check index " << endl;

  return; 
}

//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
double Electron_Level_HeI_Triplet::Ni_NeNc_LTE(double TM) const
{ 
  // HeI --> gc=2 --> ge*gc=4  
  return gw/4.0*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*TM, -1.5)
               *exp(Eion_ergs/const_kB/TM );
}

double Electron_Level_HeI_Triplet::Xi_Saha(double Xe, double Xc, double NH, double TM) const
{ return Xe*Xc*NH*Ni_NeNc_LTE(TM); }

double Electron_Level_HeI_Triplet::Ni_Saha(double Ne, double Nc, double TM) const
{ return Ne*Nc*Ni_NeNc_LTE(TM); }



//###################################################################################
//
// class: Electron_Level_HeI_Triplet_no_j
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Electron_Level_HeI_Triplet_no_j
//===================================================================================
void Electron_Level_HeI_Triplet_no_j::init(int n, int l, int mflag)
{
    nn=n; ll=l;
    mess_flag=mflag;
    gw=g_l();
    A_values.clear();
    
    //===============================================================================
    // energies
    //===============================================================================
    // quantum defect values
    if(ll<=6 && nn>10 && nn<=30)
    { 
        double dum;
        read_level_energy_QD(nameE_QD, nn, ll, 1, ll+1, dum);
        DE=dum*(2*(ll+1)+1)/(2*ll+1);
        if(l>0)
        {
            read_level_energy_QD(nameE_QD, nn, ll, 1, ll, dum);
            DE+=dum;
            read_level_energy_QD(nameE_QD, nn, ll, 1, ll-1, dum);
            DE+=dum*(2*(ll-1)+1)/(2*ll+1);
        }
        DE/=3.0;
    }
    
    // hydrogenic energies
    else if((ll>6 && nn>10) || (ll<=6 && nn>30)) DE=He1s2_ion-const_EHeI/nn/nn;
    else{ cout << " n<=10 is not allowed !" << nn << endl; exit(0); }
    //===============================================================================
    
    if(DE==-6000.0){ Dnu=Eion=Eion_ergs=nuion=0.0; cout << " BAAAADDD " << endl; }
    else
    {
        Dnu=DE*const_e/const_h;
        Eion=He1s2_ion-DE;
        Eion_ergs=Eion*const_e;
        nuion=Eion_ergs/const_h;
        
        //===========================================================================
        // all transitions from that level
        //===========================================================================
        if(HeI_Atom_read_transition_data==1)
        {
            // hydrogenic
            read_transition_inf_J_non_j(name_high_T_no_j, nn, ll, 1, -10, A_values);
        }
    }
    
    create_Transition_lookup_table();

    if(mess_flag>0)
    {
        cout << "\n %##############################################################%" << endl;
        
        display_level_information();
        cout << endl;
        display_all_downward_transitions();
        
        cout << " %##############################################################%" << endl;
    }
    
    return;
}

Electron_Level_HeI_Triplet_no_j::Electron_Level_HeI_Triplet_no_j(int n, int l, int mflag)
{ init(n, l, mflag); }

Electron_Level_HeI_Triplet_no_j::~Electron_Level_HeI_Triplet_no_j()
{ A_values.clear(); }

void Electron_Level_HeI_Triplet_no_j::create_Transition_lookup_table()
{
    ZERO_Data.np=ZERO_Data.sp=ZERO_Data.lp=ZERO_Data.jp=ZERO_Data.gwp=0;
    ZERO_Data.A21=ZERO_Data.lambda21=ZERO_Data.Dnu=0.0;
    
    Transition_Data_HeI_A_lookup_table.clear();
    vector<int> dum1;
    vector<int> dum2;
    vector<vector<int> > dum22;
    
    //=============================================================
    // get all transitions to s=0
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==0) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    int nmin=100000, nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    //=============================================================
    // get all transitions to s=1
    //=============================================================
    dum1.clear();
    for(int m=0; m<Get_n_down(); m++)
        if(A_values[m].sp==1) dum1.push_back(m);
    
    //=============================================================
    // sort according to np
    //=============================================================
    nmin=100000; nmax=0;
    for(int m=0; m<(int)dum1.size(); m++)
    {   
        if(A_values[dum1[m]].np<nmin) nmin=A_values[dum1[m]].np;
        if(A_values[dum1[m]].np>nmax) nmax=A_values[dum1[m]].np;
    }
    if(dum1.size()==0) nmax=nmin=0;
    
    dum2.clear();
    dum2.push_back(-1); // n=0..nmin transtion
    dum22.clear();
    for(int np=0; np<nmin; np++) dum22.push_back(dum2);
    
    if(nmin>=1) for(int np=nmin; np<=nmax; np++)
    {
        dum2.clear();
        
        for(int m=0; m<(int)dum1.size(); m++)
            if(A_values[dum1[m]].np==np) dum2.push_back(dum1[m]);
        
        if(dum2.size()==0) dum2.push_back(-1); 
        
        dum22.push_back(dum2);
    }       
    Transition_Data_HeI_A_lookup_table.push_back(dum22);
    
    // check the data
    if(check_Transition_Data_Tables)
    {
        Transition_Data_HeI_A T;
        for(int m=0; m<Get_n_down(); m++)
        {
            T=Get_Trans_Data(A_values[m].np, A_values[m].lp, A_values[m].sp, A_values[m].jp);
//          cout << m << " " << T.A21 << " " << A_values[m].A21 << endl;

            if(T.np==0 || T.A21-A_values[m].A21!=0.0) 
                wait_f_r(" create_Transition_lookup_table::NOOOO");
            if(T.np>nn)
            { 
                cout << " create_Transition_lookup_table::Triplet_no_j:"
                     << " there is some transtions with np>nn. Here np= " << T.np 
                     << " and nn= " << nn << endl;
                wait_f_r(); 
            }
        }
    }
    
    return;
}

void Electron_Level_HeI_Triplet_no_j::display_Trans_Data(const Transition_Data_HeI_A &v)
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet_no_j::display_Trans_Data:\n % (" 
       << nn << ", " << ll << ", 1, -10) --> (";
  // is the final Level j-resolved?
  if(v.jp!=-10) cout << v.np << ", " << v.lp << ", " << v.sp << ", " << v.jp << ") " << endl;
  else cout << v.np << ", " << v.lp << ", " << v.sp << "), gw:" << gw << endl;
  cout << " %==============================================================%" << endl;
  cout << scientific << " Dnu: " << v.Dnu*1.0e-9 << " GHz, DE: " << v.DE << " eV, lambda21: " 
       << v.lambda21 << " cm"<< endl;
  //cout << " A21: " << v.A21  << " sec^-1, f: " << v.f << endl << endl; 
  cout << scientific << " A21: " << v.A21  << " sec^-1, gw: " << v.gwp << endl << endl; 
  
  return;
}

void Electron_Level_HeI_Triplet_no_j::display_level_information()
{
  cout.precision(16);

  cout << " %==============================================================%" << endl;
  cout << " % Electron_Level_HeI_Triplet_no_j::display_level_information:\n % (n, l, s) == (" 
       << nn << ", " << ll << ", 1), gw:" << gw  
       << endl;
  cout << " %==============================================================%" << endl;
  
  cout << "\n gw: " << gw << endl;
  cout << scientific << " Dnu_1s: " << Dnu*1.0e-9 << " GHz, nuion: " 
       << nuion*1.0e-9 << " GHz" << endl;
  cout << scientific << " DE_1s: " << DE << " eV, Eion: " << Eion << " eV " << endl;
  cout << " Number of transitions: " << Get_n_down() << endl << endl;
  cout << " %==============================================================%" << endl;

  return;
}

void Electron_Level_HeI_Triplet_no_j::display_downward_transition(int i)
{
    if(i>=(int)A_values.size())
    { 
        cout << " Electron_Level_HeI_Triplet_no_j::display_downward_transition::"
             << " no transition tabulated " << i << endl;
        return;
    }
    
    display_Trans_Data(A_values[i]);

    return;
}

void Electron_Level_HeI_Triplet_no_j::display_all_downward_transitions()
{
    for(int l=0; l<(int)A_values.size(); l++) display_Trans_Data(A_values[l]);

    return;
}

//======================================================================
const Transition_Data_HeI_A& Electron_Level_HeI_Triplet_no_j::Get_Trans_Data(int np, int lp, 
                                                                             int sp, int jp) const
{
    if(sp==0 && (int)Transition_Data_HeI_A_lookup_table[0].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[0][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[0][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[0][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[0][np][m]];      
            }
        }
    }
    if(sp==1 && (int)Transition_Data_HeI_A_lookup_table[1].size()>np)
    {
        int ntrans=Transition_Data_HeI_A_lookup_table[1][np].size();
        
        for(int m=0; m<ntrans; m++) 
        {
            if(Transition_Data_HeI_A_lookup_table[1][np][m]!=-1)
            {
                if(np==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].np && 
                   sp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].sp && 
                   lp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].lp && 
                   jp==A_values[Transition_Data_HeI_A_lookup_table[1][np][m]].jp)
                    return A_values[Transition_Data_HeI_A_lookup_table[1][np][m]];      
            }
        }
    }
    
    return ZERO_Data;
}

//======================================================================
double Electron_Level_HeI_Triplet_no_j::Get_A21(int np, int lp, int sp, int jp) const
{
    if(np<=(int)HeI_Atom_njresolved) 
        for(int i=0; i<(int)A_values.size(); i++) 
        {
            if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
                return A_values[i].A21;
            //else return 0.0;
        }
    else 
        for(int i=0; i<(int)A_values.size(); i++) 
        {
            if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==-10) 
                return A_values[i].A21;
            //else return 0.0;
        }
    
  return 0.0;
}

/*
double Electron_Level_HeI_Triplet_no_j::Get_f(int np, int lp, int sp, int jp) const
{
  if(np<=(int)HeI_Atom_njresolved) 
    for(int i=0; i<(int)A_values.size(); i++) 
      {
    if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
        return A_values[i].f;
    else return 0.0;
      }
  else 
    for(int i=0; i<A_values.size(); i++) 
      {
    if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==-10) 
        return A_values[i].f;
    else return 0.0;
      }

  return 0.0;
}
*/

double Electron_Level_HeI_Triplet_no_j::Get_nu21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].Dnu;
        //else return 0.0;
    }
    
    return 0.0;
}

double Electron_Level_HeI_Triplet_no_j::Get_lambda21(int np, int lp, int sp, int jp) const
{
    for(int i=0; i<(int)A_values.size(); i++) 
    {
        if(A_values[i].np==np && A_values[i].lp==lp && A_values[i].sp==sp && A_values[i].jp==jp) 
            return A_values[i].lambda21;
        //else return 0.0;
    }
    
    return 0.0;
}

void Electron_Level_HeI_Triplet_no_j::Set_xxx_Transition_Data(int m, int d, double v) 
{ 
  if(d<=5)
    {
      A_values[m].xxx[d]=v; 
    }
  else cout << " Electron_Level_HeI_Triplet_no_j::Set_xxx_Transition_Data: check index " << endl;

  return; 
}

//==============================================================================
// Saha-relations with continuum
//==============================================================================
// f(T) == (Ni/[Ne Nc])_LTE
double Electron_Level_HeI_Triplet_no_j::Ni_NeNc_LTE(double TM) const
{ 
  // HeI --> gc=2 --> ge*gc=4  
  return gw/4.0*pow(const_lambdac, 3)*pow(2.0*PI*const_kb_mec2*TM, -1.5)
               *exp(Eion_ergs/const_kB/TM );
}

double Electron_Level_HeI_Triplet_no_j::Xi_Saha(double Xe, double Xc, double NH, double TM) const
{ return Xe*Xc*NH*Ni_NeNc_LTE(TM); }

double Electron_Level_HeI_Triplet_no_j::Ni_Saha(double Ne, double Nc, double TM) const
{ return Ne*Nc*Ni_NeNc_LTE(TM); }

//###################################################################################
//
// class: Atomic_Shell_HeI_Singlet
//
//###################################################################################

//===================================================================================
//Konstructors and Destructors for class: Atomic_Shell_HeI_Singlet
//===================================================================================
void Atomic_Shell_HeI_Singlet::init(int n, int mflag)
{
    nn=n;
    mess_flag=mflag;
    create_Electron_Levels();
}

Atomic_Shell_HeI_Singlet::Atomic_Shell_HeI_Singlet(int n, int mflag)
{ init(n, mflag); }
    
Atomic_Shell_HeI_Singlet::~Atomic_Shell_HeI_Singlet()
{ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell_HeI_Singlet::create_Electron_Levels()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell_HeI_Singlet::create_Electron_Levels: filling shell: " 
             << nn <<"\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();
    
    Electron_Level_HeI_Singlet v;
    // fill with empty electron-states
    for(int l=0; l<nn; l++) Angular_Momentum_Level.push_back(v);
    // now initialize each state
    for(int l=0; l<nn; l++) Angular_Momentum_Level[l].init(nn, l, mess_flag);
    
    return;
}

void Atomic_Shell_HeI_Singlet::display_general_data_of_level(int i)
{
    if(i>=(int)(nn) || Angular_Momentum_Level.size()==0)
    {
    cout << " This level does not exist inside shell " << nn << endl;
    return;
    }

    cout << " %==============================================================%" << endl;
    cout << " % Atomic_Shell_HeI_Singlet::display_general_data_of_level:\n %\n" << endl; 
  
    Angular_Momentum_Level[i].display_level_information();

    return;
}

//###################################################################################
//
//class: Atomic_Shell_HeI_Triplet
//
//###################################################################################
//===================================================================================
// Konstructors and Destructors for class: Atomic_Shell_HeI_Triplet
//===================================================================================
void Atomic_Shell_HeI_Triplet::init(int n, int mflag)
{
    nn=n;
    mess_flag=mflag;
    create_Electron_Levels();
}

Atomic_Shell_HeI_Triplet::Atomic_Shell_HeI_Triplet(int n, int mflag)
{ init(n, mflag); }
    
Atomic_Shell_HeI_Triplet::~Atomic_Shell_HeI_Triplet()
{ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell_HeI_Triplet::create_Electron_Levels()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell_HeI_Triplet::create_Electron_Levels: filling shell: " 
             << nn <<"\n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();

    Electron_Level_HeI_Triplet v;
    vector<Electron_Level_HeI_Triplet> v3;
    v3.push_back(v);
    v3.push_back(v);
    v3.push_back(v);

    // fill with empty electron-states
    for(int l=0; l<nn; l++) Angular_Momentum_Level.push_back(v3);
    // now initialize each state
    // n3S -state
    Angular_Momentum_Level[0][2].init(nn, 0, 1, mess_flag);
    // n3X -states
    for(int l=1; l<nn; l++)
      {
          Angular_Momentum_Level[l][0].init(nn, l, l-1, mess_flag);
          Angular_Momentum_Level[l][1].init(nn, l, l, mess_flag);
          Angular_Momentum_Level[l][2].init(nn, l, l+1, mess_flag);
      }

    return;
}

void Atomic_Shell_HeI_Triplet::display_general_data_of_level(int i)
{
    if(i>=(int)(nn) || Angular_Momentum_Level.size()==0)
    {
        cout << " This level does not exist inside shell " << nn << endl;
        return;
    }
    
    cout << " %==============================================================%" << endl;
    cout << " % Atomic_Shell_HeI_Triplet::display_general_data_of_level:\n %\n" << endl; 
    
    Angular_Momentum_Level[i][2].display_level_information();
    if(i!=0)  
    {
        Angular_Momentum_Level[i][1].display_level_information();
        Angular_Momentum_Level[i][0].display_level_information();
    }
    return;
}

//###################################################################################
//
// class: Atomic_Shell_HeI_Triplet_no_j
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for class: Atomic_Shell_HeI_Triplet_no_j
//===================================================================================
void Atomic_Shell_HeI_Triplet_no_j::init(int n, int mflag)
{
    nn=n;
    mess_flag=mflag;
    create_Electron_Levels();
}

Atomic_Shell_HeI_Triplet_no_j::Atomic_Shell_HeI_Triplet_no_j(int n, int mflag)
{ init(n, mflag); }
    
Atomic_Shell_HeI_Triplet_no_j::~Atomic_Shell_HeI_Triplet_no_j()
{ Angular_Momentum_Level.clear(); }

//-----------------------------------------------------------------------------------
void Atomic_Shell_HeI_Triplet_no_j::create_Electron_Levels()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atomic_Shell_HeI_Triplet_no_j::create_Electron_Levels: filling shell: " 
             << nn <<"\n % " << endl;
        cout << " % The levels will not be j-resolved \n % " << endl;
        cout << " %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    // free memory (if necessary)
    Angular_Momentum_Level.clear();

    Electron_Level_HeI_Triplet_no_j v;
    // fill with empty electron-states
    for(int l=0; l<nn; l++) Angular_Momentum_Level.push_back(v);
    // now initialize each state
    for(int l=0; l<nn; l++) Angular_Momentum_Level[l].init(nn, l, mess_flag);

    return;
}

void Atomic_Shell_HeI_Triplet_no_j::display_general_data_of_level(int i)
{
    if(i>=(int)(nn) || Angular_Momentum_Level.size()==0)
    {
        cout << " This level does not exist inside shell " << nn << endl;
        return;
    }
    
    cout << " %==============================================================%" << endl;
    cout << " % Atomic_Shell_HeI_Triplet_no_j::display_general_data_of_level:\n %\n" << endl; 
    
    Angular_Momentum_Level[i].display_level_information();
    
    return;
}

//###################################################################################
//
// class: Atom_HeI_Singlet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Atom_HeI_Singlet
//===================================================================================
void Atom_HeI_Singlet::init(int nS, int mflag)
{
    nShells=nS;
    mess_flag=mflag;

    fill_Level_Map();

    create_Shells();
    if(mess_flag>2) display_Level_Map();
}

Atom_HeI_Singlet:: Atom_HeI_Singlet(){}

Atom_HeI_Singlet:: Atom_HeI_Singlet(int nS, int mflag){ init(nS, mflag); }

Atom_HeI_Singlet::~Atom_HeI_Singlet()
{ 
  Level_Map.clear(); 
  Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom_HeI_Singlet::create_Shells()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom_HeI_Singlet::create_Shells:creating all the shells up to n=" 
             << nShells << endl;
        cout << " %\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();

    Atomic_Shell_HeI_Singlet v;
    // fill with empty shells
    for(int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(int n=1; n<=nShells; n++) Shell[n].init(n, m); 

    return;
}

//------------------------------------------------------------------------------------------------------
void Atom_HeI_Singlet::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 
    
    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(int n=1; n<=nShells; n++) 
        for(int l=0; l<n; l++)
        {
            v.n=n; v.l=l;         
            Level_Map.push_back(v);
        }
    
    return;
}

void Atom_HeI_Singlet::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom_HeI_Singlet::display_Level_Map:" << endl; 
    cout << " % Format: i [check of inversion] --> (n, l, s, j)" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(int i=0; i<(int)Level_Map.size(); i++)
    {
        if(Get_Level_index(Level_Map[i].n, 0)==(int)i) cout << " Shell " << Level_Map[i].n << endl;
        
        cout << " " << i << " [" << Get_Level_index(Level_Map[i].n, Level_Map[i].l) << "] --> (" 
             << Level_Map[i].n << ", " << Level_Map[i].l << ", 0, " 
             << Level_Map[i].l <<")" << endl; 
    }
    cout << endl;
    
    return;
}

//================================================================================
// to check level access
//================================================================================
void Atom_HeI_Singlet::check_Level(int n, int l, string mess) const
{
    if(n<1 || n>nShells || l>=n || l<0)
    {
        cout << " Atom_HeI_Singlet::" << mess << ": you are trying to access a non-existing level: " << n 
             << ", " << l << ", 0, " << l << endl;
        exit(0);
    }
    return;
}

void Atom_HeI_Singlet::check_Level(int i, string mess) const
{
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom_HeI_Singlet::" << mess << ": you are trying to access a non-existing level: " << i 
             << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    return;
}

const Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int i) const
{ 
    check_Level(i, "Level");

    //cout << " Atom_HeI_Singlet: " << i << " " << Level_Map[i].n << " " << Level_Map[i].l << endl;

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}
 
const Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int n, int l) const
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}
 
Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int i) 
{ 
   check_Level(i, "Level");

   //cout << " Atom_HeI_Singlet: " << i << " " << Level_Map[i].n << " " << Level_Map[i].l << endl;

   return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}
 
Electron_Level_HeI_Singlet& Atom_HeI_Singlet::Level(int n, int l)
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}
 
//###################################################################################
//
// class: Atom_HeI_Triplet
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Atom_HeI_Triplet
//===================================================================================
void Atom_HeI_Triplet::init(int nS, int mflag)
{
    nShells=nS;
    mess_flag=mflag;

    fill_Level_Map();
    create_Shells();

    if(mess_flag>2) display_Level_Map();
}

Atom_HeI_Triplet:: Atom_HeI_Triplet(){}

Atom_HeI_Triplet:: Atom_HeI_Triplet(int nS, int mflag){ init(nS, mflag); }

Atom_HeI_Triplet::~Atom_HeI_Triplet()
{ 
  Level_Map.clear(); 
  Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom_HeI_Triplet::create_Shells()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom_HeI_Triplet::create_Shells:creating all the shells up to n=" 
             << nShells << endl;
        cout << " %\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();
    
    Atomic_Shell_HeI_Triplet v;
    // fill with empty shells
    for(int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(int n=2; n<=nShells; n++) Shell[n].init(n, m); 
    
    return;
}

int Atom_HeI_Triplet::Get_Level_index(int n, int l, int j) const 
{ 
    if(n<2 || n>nShells || l>=n || l<0 || l+1<j || l>j+1 || (l==0 && j!=1)) 
        cout << " Atom_HeI_Triplet::Get_Level_index: no level with " 
             << n << " " << l << " " << j << " !!! " << endl;

  return Get_number_of_Levels_until(n-1)+3*(l+1)-2-(l+1-j)-1; 
}

//------------------------------------------------------------------------------------------------------
void Atom_HeI_Triplet::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 
    
    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(int n=2; n<=nShells; n++) 
    {
        v.n=n; v.l=0; v.j=1;      
        Level_Map.push_back(v);
        
        for(int l=1; l<n; l++)
        {
            v.n=n; v.l=l; v.j=l-1;        
            Level_Map.push_back(v);
            v.n=n; v.l=l; v.j=l;      
            Level_Map.push_back(v);
            v.n=n; v.l=l; v.j=l+1;        
            Level_Map.push_back(v);
        }
    }
    
    return;
}

void Atom_HeI_Triplet::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom_HeI_Triplet::display_Level_Map:" << endl; 
    cout << " % Format: i --> (n, l, s, j)" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(int i=0; i<(int)Level_Map.size(); i++)
    {
        if(Get_Level_index(Level_Map[i].n, 0, 1)==(int)i) cout << " Shell " << Level_Map[i].n << endl;
        
        cout << " " << i << " [" << Get_Level_index(Level_Map[i].n, Level_Map[i].l, Level_Map[i].j) << "] --> (" 
             << Level_Map[i].n << ", " << Level_Map[i].l << ", 1, " 
             << Level_Map[i].j <<")" << endl; 
    }
    
    cout << endl;
    
    return;
}

//================================================================================
// to check level access
//================================================================================
void Atom_HeI_Triplet::check_Level(int n, int l, int j, string mess) const
{
    if(n<2 || n>nShells || l>=n || l<0 || l+1<j || l>j+1 || (l==0 && j!=1))
    {
        cout << " Atom_HeI_Triplet::" << mess 
             << ": you are trying to access a non-existing level: " << n 
             << ", " << l << ", 1, " << j << endl;
        exit(0);
    }
    return;
}

void Atom_HeI_Triplet::check_Level(int i, string mess) const
{
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom_HeI_Triplet::" << mess 
             << ": you are trying to access a non-existing level: " << i 
             << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    return;
}

const Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int i) const
{ 
    check_Level(i, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l][Level_Map[i].j-(Level_Map[i].l-1)]; 
}
 
const Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int n, int l, int j) const
{ 
    check_Level(n, l, j, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[n].Angular_Momentum_Level[l][j-(l-1)]; 
}
 
Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int i) 
{ 
    check_Level(i, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l][Level_Map[i].j-(Level_Map[i].l-1)]; 
}
 
Electron_Level_HeI_Triplet& Atom_HeI_Triplet::Level(int n, int l, int j)
{ 
    check_Level(n, l, j, "Level");

    // the access of data is such that for l==0 --> j-(l-1)=2
    return Shell[n].Angular_Momentum_Level[l][j-(l-1)]; 
}

//###################################################################################
//
// class: Atom_HeI_Triplet_no_j
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for Atom_HeI_Triplet_no_j
//===================================================================================
void Atom_HeI_Triplet_no_j::init(int nS, int mflag)
{
    nShells=nS;
    mess_flag=mflag;

    fill_Level_Map();
    create_Shells();

    if(mess_flag>2) display_Level_Map();
}

Atom_HeI_Triplet_no_j:: Atom_HeI_Triplet_no_j(){}

Atom_HeI_Triplet_no_j:: Atom_HeI_Triplet_no_j(int nS, int mflag){ init(nS, mflag); }

Atom_HeI_Triplet_no_j::~Atom_HeI_Triplet_no_j()
{ 
  Level_Map.clear(); 
  Shell.clear();
}

//-----------------------------------------------------------------------------------
void Atom_HeI_Triplet_no_j::create_Shells()
{
    if(mess_flag>-1)
    {
        cout << "\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n %" << endl;
        cout << " % Atom_HeI_Triplet_no_j::create_Shells:creating all the shells from " 
             << HeI_Atom_njresolved+1 << " up to n=" << nShells << endl;
        cout << " %\n %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%\n" << endl;
    }
    
    int m=mess_flag-1;
    if(mess_flag>=1) m++;
    if(mess_flag>=2) m++;
    
    // free memory (if necessary)
    Shell.clear();
    
    Atomic_Shell_HeI_Triplet_no_j v;
    // fill with empty shells
    for(int n=0; n<=nShells; n++) Shell.push_back(v);
    // create each shells
    for(int n=HeI_Atom_njresolved+1; n<=nShells; n++) Shell[n].init(n, m); 
    
    return;
}

int Atom_HeI_Triplet_no_j::Get_number_of_Levels_until(int nmax) const
{ return nmax*(nmax+1)/2 - HeI_Atom_njresolved*(HeI_Atom_njresolved+1)/2; }

int Atom_HeI_Triplet_no_j::Get_Level_index(int n, int l) const 
{ 
    if(n>nShells || n<=HeI_Atom_njresolved || l>=n || l<0) 
        cout << " Atom_HeI_Triplet_no_j::Get_Level_index: no level with " 
             << n << " " << l << " !!! " << endl;
    
    return Get_number_of_Levels_until(n-1)+l; 
}

//------------------------------------------------------------------------------------------------------
void Atom_HeI_Triplet_no_j::fill_Level_Map()
{
    // free memory (if necessary)
    Level_Map.clear(); 

    // this function creates the map of indicies i-> (n,l)
    Level_nl v;
    
    for(int n=HeI_Atom_njresolved+1; n<=nShells; n++) 
      for(int l=0; l<n; l++)
    {
      v.n=n; v.l=l;     
      Level_Map.push_back(v);
    }
   
    return;
}

void Atom_HeI_Triplet_no_j::display_Level_Map()
{
    cout << "\n %==============================================================%" << endl;
    cout << " % Atom_HeI_Triplet_no_j::display_Level_Map:" << endl; 
    cout << " % Format: i --> (n, l, s, -10)" << endl; 
    cout << " %==============================================================%" << endl;
    
    for(int i=0; i<(int)Level_Map.size(); i++)
    {
        if(Get_Level_index(Level_Map[i].n, 0)==(int)i) cout << " Shell " << Level_Map[i].n << endl;
        
        cout << " " << i << " [" << Get_Level_index(Level_Map[i].n, Level_Map[i].l) << "] --> (" 
             << Level_Map[i].n << ", " << Level_Map[i].l << ", 1, -10)" << endl; 
    }
    
    cout << endl;
    
    return;
}

//================================================================================
// to check level access
//================================================================================
void Atom_HeI_Triplet_no_j::check_Level(int n, int l, string mess) const
{
    if(n>nShells || n<=HeI_Atom_njresolved || l>=n || l<0)
    {
        cout << " Atom_HeI_Triplet_no_j::" << mess 
             << ": you are trying to access a non-existing level: " << n 
             << ", " << l << ", 1 " << endl;
        exit(0);
    }
    return;
}

void Atom_HeI_Triplet_no_j::check_Level(int i, string mess) const
{
    if(i> Get_total_number_of_Levels())
    {
        cout << " Atom_HeI_Triplet_no_j::" << mess 
             << ": you are trying to access a non-existing level: " << i 
             << " total number of levels: " << Get_total_number_of_Levels() << endl;
        exit(0);
    }
    return;
}

const Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int i) const
{ 
    check_Level(i, "Level");

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l]; 
}
 
const Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int n, int l) const
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}
 
Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int i) 
{ 
    check_Level(i, "Level");

    return Shell[Level_Map[i].n].Angular_Momentum_Level[Level_Map[i].l];
}
 
Electron_Level_HeI_Triplet_no_j& Atom_HeI_Triplet_no_j::Level(int n, int l)
{ 
    check_Level(n, l, "Level");

    return Shell[n].Angular_Momentum_Level[l]; 
}

//###################################################################################
//
// class: Gas_of_HeI_Atoms
//
//###################################################################################

//===================================================================================
// Konstructors and Destructors for class: Gas_of_HeI_Atoms
//===================================================================================
void Gas_of_HeI_Atoms::init(int nS, int njres, int mflag)
{ 
    HeI_Atom_njresolved=(int)min(njres, nS);
    // do not read transition data
    if(nS<0){ HeI_Atom_read_transition_data=0; nS=-nS; }
    
    Sing.init(nS, mflag);
    Trip.init((int)min(nS, HeI_Atom_njresolved), mflag);
    if(nS>(int)HeI_Atom_njresolved) Trip_no_j.init(nS, mflag);
    
    indexT=Sing.Get_total_number_of_Levels();
    nl=indexT_no_j=indexT+Trip.Get_total_number_of_Levels();
    if(nS>(int)HeI_Atom_njresolved) nl+=Trip_no_j.Get_total_number_of_Levels(); 
    
    if(HeI_Atom_read_transition_data!=0) 
    {
        //============================================================
        // check the consistency of transition data (17.05.2009)
        //============================================================
        check_transition_data();        
        
        //===========================================================
        // create the Voigt profiles
        //===========================================================
        double A21, f=0.0, nu21, lam21, Gamma;
       
        //===========================================================
        // Transition-Data for n^1 P_1 - 1^1 S_0
        //===========================================================
        for(int n=2; n<=(int)min(nS, 10); n++)
        {
            A21=Sing.Level(n, 1).Get_A21(1, 0, 0, 0);
            nu21=Sing.Level(n, 1).Get_nu21(1, 0, 0, 0);
            lam21=Sing.Level(n, 1).Get_lambda21(1, 0, 0, 0);
            
            Gamma=0.0;
            for(int m=0; m<(int)Sing.Level(n, 1).Get_n_down(); m++) 
                Gamma+=Sing.Level(n, 1).Get_Trans_Data(m).A21;
            
            if(mflag>=1) cout << " Initializing profile for " << n 
                              << "^1 P_1 - 1^1 S_0 transition:" 
                              << nu21 << " " << lam21 << " " << A21 << " " 
                              << Gamma << " " << f << endl;

            phi_HeI_nP_S[n].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 4);
        }
        
        //===========================================================
        // Transition-Data for n^3 P_1 - 1^1 S_0
        //===========================================================
        for(int n=2; n<=(int)min(nS, _HeI_add_Intercomb); n++)
        {
            A21=Trip.Level(n, 1, 1).Get_A21(1, 0, 0, 0);
            nu21=Trip.Level(n, 1, 1).Get_nu21(1, 0, 0, 0);
            lam21=Trip.Level(n, 1, 1).Get_lambda21(1, 0, 0, 0);
            
            Gamma=0.0;
            for(int m=0; m<(int)Trip.Level(n, 1, 1).Get_n_down(); m++) 
                Gamma+=Trip.Level(n, 1, 1).Get_Trans_Data(m).A21;
            
            if(mflag>=1) cout << " Initializing profile for " << n 
                              << "^3 P_1 - 1^1 S_0 transition:" 
                              << nu21 << " " << lam21 << " " << A21 << " " 
                              << Gamma << " " << f << endl;

            phi_HeI_nP_T[n].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 4);
        }

        //===========================================================
        // Transition-Data for n^1 D_1 - 1^1 S_0
        //===========================================================
        for(int n=3; n<=(int)min(nS, _HeI_add_Quad); n++)
        {
            A21=Sing.Level(n, 2).Get_A21(1, 0, 0, 0);
            nu21=Sing.Level(n, 2).Get_nu21(1, 0, 0, 0);
            lam21=Sing.Level(n, 2).Get_lambda21(1, 0, 0, 0);
            
            Gamma=0.0;
            for(int m=0; m<(int)Sing.Level(n, 2).Get_n_down(); m++) 
                Gamma+=Sing.Level(n, 2).Get_Trans_Data(m).A21;
            
            if(mflag>=1) cout << " Initializing profile for " << n 
                              << "^1 D_1 - 1^1 S_0 transition:" 
                              << nu21 << " " << lam21 << " " << A21 << " " 
                              << Gamma << " " << f << endl;
            
            phi_HeI_nD_S[n].Set_atomic_data_Gamma(nu21, lam21, A21, f, Gamma, 4);
        }
    }
    
    //========================================
    // create HI Ly-c cross section
    //========================================
    HILyc.init(1, 0, 1, 1, 0);

    return;
}

Gas_of_HeI_Atoms::Gas_of_HeI_Atoms(){}

Gas_of_HeI_Atoms::Gas_of_HeI_Atoms(int nS, int njres, int mflag)
{ init(nS, njres, mflag); }
    
Gas_of_HeI_Atoms::~Gas_of_HeI_Atoms(){}

//===================================================================================
// Access to Voigt-profiles
//===================================================================================
const Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_S_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), 10)) return phi_HeI_nP_S[n]; 
    
    return phi_HeI_nP_S[0]; 
}

Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_S_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), 10)) return phi_HeI_nP_S[n]; 
    
    return phi_HeI_nP_S[0]; 
}

//===================================================================================
const Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_T_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), _HeI_add_Intercomb)) return phi_HeI_nP_T[n]; 
    
    return phi_HeI_nP_T[0]; 
}

Voigtprofile_Dawson& Gas_of_HeI_Atoms::nP_T_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), _HeI_add_Intercomb)) return phi_HeI_nP_T[n]; 
    
    return phi_HeI_nP_T[0]; 
}

//===================================================================================
const Voigtprofile_Dawson& Gas_of_HeI_Atoms::nD_S_profile(int n) const 
{ 
    if(n>1 && n<=(int)min(Get_nShells(), _HeI_add_Quad)) return phi_HeI_nD_S[n]; 
    
    return phi_HeI_nD_S[0]; 
}

Voigtprofile_Dawson& Gas_of_HeI_Atoms::nD_S_profile(int n)
{ 
    if(n>1 && n<=(int)min(Get_nShells(), _HeI_add_Quad)) return phi_HeI_nD_S[n]; 
    
    return phi_HeI_nD_S[0]; 
}

//===================================================================================
void Gas_of_HeI_Atoms:: check_transition_data()
{
    double nucu, nucl, Dnu, val;
    int ip, count=0;
    
    for(int k=0; k<(int)nl; k++)
    {
        //cout << "\n Level " << k << " == (" << Get_n(k) << ", " <<  Get_l(k) << ", " << Get_S(k) << ", " <<  Get_J(k) << ")" << endl;
        //cout << " Total number of down transitions = " << Get_n_down(k) << endl;
        
        nucu=Get_nu_ion(k);
        for(int m=0; m<(int)Get_n_down(k); m++)
        {
            ip=Get_Level_index(Get_Trans_Data(k, m).np, Get_Trans_Data(k, m).lp, Get_Trans_Data(k, m).sp, Get_Trans_Data(k, m).jp);
            
            nucl=Get_nu_ion(ip);
            Dnu=nucl-nucu;
            
            val=Get_Trans_Data(k, m).Dnu;
            
            if(fabs(Dnu/val-1.0)>=1.0e-8) 
            {
                count++;
                cout << "\n Level " << k << " == (" << Get_n(k) << ", " <<  Get_l(k) << ", " << Get_S(k) << ", " <<  Get_J(k) << ")" << endl;
                cout << " There is an inconsistency of " << fabs(Dnu/val-1.0) << " in the transition frequency " << val << " Hz to level " 
                << ip << " == (" << Get_Trans_Data(k, m).np << ", " << Get_Trans_Data(k, m).lp << ", " << Get_Trans_Data(k, m).sp << ", " << Get_Trans_Data(k, m).jp << ")" << endl; 
            }
        }
    }
    
    if(count!=0)
    {
        cout << " There were " << count << " inconsistencies in the transition frequencies. Maybe you want to recompute your Helium model... " << endl;
        wait_f_r();
    }
    
    return;
}

//===================================================================================
int Gas_of_HeI_Atoms::Get_njres(){ return HeI_Atom_njresolved; }

const Transition_Data_HeI_A& Gas_of_HeI_Atoms::Get_Trans_Data(int i, int np, int lp, 
                                                              int sp, int jp) const
{
    int n=Get_n(i);
    int l=Get_l(i);
    int j=Get_J(i);
    int s=Get_S(i);
    
    if(s==0) return Sing.Level(n, l).Get_Trans_Data(np, lp, sp, jp);
    else if(s==1 && n<=(int)HeI_Atom_njresolved) 
        return Trip.Level(n, l, j).Get_Trans_Data(np, lp, sp, jp);
    else if(s==1 && n>(int)HeI_Atom_njresolved) 
        return Trip_no_j.Level(n, l).Get_Trans_Data(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    
    return Sing.Level(i).Get_Trans_Data(0, lp, sp, jp);
}

double Gas_of_HeI_Atoms:: Get_A(int n, int l, int s, int j, int np, int lp, int sp, int jp) const
{
  if(s==0) return Sing.Level(n, l).Get_A21(np, lp, sp, jp);
  else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Level(n, l, j).Get_A21(np, lp, sp, jp);
  else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Level(n, l).Get_A21(np, lp, sp, jp);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
    
  return 0.0;
}

/*
double Gas_of_HeI_Atoms:: Get_f(int n, int l, int s, int j, 
                                int np, int lp, int sp, int jp) const
{
    if(s==0) return Sing.Level(n, l).Get_f(np, lp, sp, jp);
    else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Level(n, l, j).Get_f(np, lp, sp, jp);
    else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Level(n, l).Get_f(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    return 0.0;
}
*/

double Gas_of_HeI_Atoms:: Get_nu21(int n, int l, int s, int j, 
                                   int np, int lp, int sp, int jp) const
{
    if(s==0) return Sing.Level(n, l).Get_nu21(np, lp, sp, jp);
    else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Level(n, l, j).Get_nu21(np, lp, sp, jp);
    else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Level(n, l).Get_nu21(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    return 0.0;
}

double Gas_of_HeI_Atoms:: Get_lambda21(int n, int l, int s, int j, 
                                       int np, int lp, int sp, int jp) const
{
    if(s==0) return Sing.Level(n, l).Get_lambda21(np, lp, sp, jp);
    else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Level(n, l, j).Get_lambda21(np, lp, sp, jp);
    else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Level(n, l).Get_lambda21(np, lp, sp, jp);
    
    cout << " Gas_of_HeI_Atoms:: fail " << endl;
    return 0.0;
}

double Gas_of_HeI_Atoms::Xi(int n, int l, int s, int j) const
{
  if(s==0) return Sing.Level(n, l).Get_Xi();
  else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Level(n, l, j).Get_Xi();
  else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Level(n, l).Get_Xi();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Ric(int n, int l, int s, int j) const
{
  if(s==0) return Sing.Level(n, l).Get_Ric();
  else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Level(n, l, j).Get_Ric();
  else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Level(n, l).Get_Ric();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

void Gas_of_HeI_Atoms::Set_Xi(int n, int l, int s, int j, double X)
{
  if(s==0) Sing.Level(n, l).Set_Xi(X);
  else if(s==1 && n<=HeI_Atom_njresolved) Trip.Level(n, l, j).Set_Xi(X);
  else if(s==1 && n>HeI_Atom_njresolved) Trip_no_j.Level(n, l).Set_Xi(X);
  return;
}

void Gas_of_HeI_Atoms::Set_Ric(int n, int l, int s, int j, double X)
{
  if(s==0) Sing.Level(n, l).Set_Ric(X);
  else if(s==1 && n<=HeI_Atom_njresolved) Trip.Level(n, l, j).Set_Ric(X);
  else if(s==1 && n>HeI_Atom_njresolved) Trip_no_j.Level(n, l).Set_Ric(X);
  return;
}

//=================================================================================================
// show level info
//=================================================================================================
void Gas_of_HeI_Atoms::Show_NLSJ(int i) const
{
    cout << " Gas_of_HeI_Atoms::Show_NLSJ: (N, L, S, J) == (" 
         << Get_n(i) << ", " << Get_l(i) << ", " 
         << Get_S(i) << ", " << Get_J(i) << ") " 
         << endl;
    
    return;
}

//=================================================================================================
// versions that follow Singlet & Triplet sub-sequently
//=================================================================================================
int Gas_of_HeI_Atoms::Get_Level_index(int n, int l, int s, int j) const
{
  if(s==0) return Sing.Get_Level_index(n, l);
  else if(s==1 && n<=HeI_Atom_njresolved) return Trip.Get_Level_index(n, l, j)+indexT;
  else if(s==1 && n>HeI_Atom_njresolved) return Trip_no_j.Get_Level_index(n, l)+indexT_no_j;

  cout << " Get_as_of_HeI_Atoms::Get_Level_index: No level with s= " << s << endl; 
  return 0; 
}

//===================================================================================
int Gas_of_HeI_Atoms::Get_n(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_n();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_n();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_n();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

int Gas_of_HeI_Atoms::Get_l(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_l();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_l();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_l();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

int Gas_of_HeI_Atoms::Get_S(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_S();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_S();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_S();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

int Gas_of_HeI_Atoms::Get_J(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_J();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_J();
  else if(i<nl) return -10;

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

double Gas_of_HeI_Atoms::Get_gw(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_gw();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_gw();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_gw();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Get_nu_ion(int i) const
{
    if(i<indexT) return Sing.Level(i).Get_nu_ion();
    else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_nu_ion();
    else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_nu_ion();
    
    cout << " Gas_of_HeI_Atoms:: Get_nu_ion: fail " << endl;
    return 0.0;
}
    
//===================================================================================
int Gas_of_HeI_Atoms::Get_n_down() const
{
  int s=0;
  for(int i=0; i<indexT; i++) s+=Sing.Level(i).Get_n_down();
  for(int i=indexT; i<indexT_no_j; i++) s+=Trip.Level(i-indexT).Get_n_down();
  for(int i=indexT_no_j; i<nl; i++) s+=Trip_no_j.Level(i-indexT_no_j).Get_n_down();

  return s;
}

int Gas_of_HeI_Atoms::Get_n_down(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_n_down();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_n_down();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_n_down();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 1;
}

const Transition_Data_HeI_A& Gas_of_HeI_Atoms::Get_Trans_Data(int i, int m) const
{
  if(i<indexT) return Sing.Level(i).Get_Trans_Data(m);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_Trans_Data(m);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_Trans_Data(m);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return Sing.Level(0).Get_Trans_Data(0);
}

//===================================================================================
double Gas_of_HeI_Atoms::Xi(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_Xi();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_Xi();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_Xi();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Ric(int i) const
{
  if(i<indexT) return Sing.Level(i).Get_Ric();
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Get_Ric();
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Get_Ric();

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

void Gas_of_HeI_Atoms::Set_Xi(int i, double X)
{
  if(i<indexT) Sing.Level(i).Set_Xi(X);
  else if(i<indexT_no_j) Trip.Level(i-indexT).Set_Xi(X);
  else if(i<nl) Trip_no_j.Level(i-indexT_no_j).Set_Xi(X);

  return;
}

void Gas_of_HeI_Atoms::Set_Ric(int i, double X)
{
  if(i<indexT) Sing.Level(i).Set_Ric(X);
  else if(i<indexT_no_j) Trip.Level(i-indexT).Set_Ric(X);
  else if(i<nl) Trip_no_j.Level(i-indexT_no_j).Set_Ric(X);

  return;
}

double Gas_of_HeI_Atoms::X_tot() const
{
    double r=0.0;   
    // smallest terms first!
    for(int i=nl-1; i>=indexT_no_j; i--) r+=Trip_no_j.Level(i-indexT_no_j).Get_Xi();
    for(int i=indexT_no_j-1; i>=indexT; i--) r+=Trip.Level(i-indexT).Get_Xi();
    // carefull here! i=-1 is non-sense; need signed int here!
    for(int i=indexT-1; i>=0; i--) r+=Sing.Level(i).Get_Xi();
    return r; 
}

void Gas_of_HeI_Atoms::Set_xxx_Transition_Data(int i, int m, int d, double v)
{
  if(i<indexT) Sing.Level(i).Set_xxx_Transition_Data(m, d, v);
  else if(i<indexT_no_j) Trip.Level(i-indexT).Set_xxx_Transition_Data(m, d, v);
  else if(i<nl) Trip_no_j.Level(i-indexT_no_j).Set_xxx_Transition_Data(m, d, v);
  return;
}

//===================================================================================
// Saha-relations with continuum
//===================================================================================
double Gas_of_HeI_Atoms::Ni_NeNc_LTE(int i, double TM) const
{
  if(i<indexT) return Sing.Level(i).Ni_NeNc_LTE(TM);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Ni_NeNc_LTE(TM);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Ni_NeNc_LTE(TM);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Xi_Saha(int i, double Xe, double Xc, double NH, double TM)  const
{
  if(i<indexT) return Sing.Level(i).Xi_Saha(Xe, Xc, NH, TM);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Xi_Saha(Xe, Xc, NH, TM);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Xi_Saha(Xe, Xc, NH, TM);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}

double Gas_of_HeI_Atoms::Ni_Saha(int i, double Ne, double Nc, double TM)  const
{
  if(i<indexT) return Sing.Level(i).Ni_Saha(Ne, Nc, TM);
  else if(i<indexT_no_j) return Trip.Level(i-indexT).Ni_Saha(Ne, Nc, TM);
  else if(i<nl) return Trip_no_j.Level(i-indexT_no_j).Ni_Saha(Ne, Nc, TM);

  cout << " Gas_of_HeI_Atoms:: fail " << endl;
  return 0.0;
}
 
