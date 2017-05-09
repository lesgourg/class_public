//========================================================================
// data for helium atom
//========================================================================

Gas_of_HeI_Atoms *HeIAtoms;
#define HeI_Atoms (*HeIAtoms)

//========================================================================
// switch quadrupole lines on/off (added 29.12.2010)
//========================================================================
bool set_HeI_Quadrupole_lines=1;

//========================================================================
//
// setting for helium effective rates
//
//========================================================================
const int nS_effective_HeI=30;                     // this is always fixed

string effective_rate_path_HeI_arr[7]=
{
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.res_2/",
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.2s2p3p_ST/",
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.res_5/",
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.res_10/",
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.2s2p3p_ST_Q/",
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.res_5_Q/",
    Rec_database_path+"/Effective_Rates.HeI/Effective_Rate_Tables.HeI.res_10_Q/"
};

string effective_rate_path_HeI;

//========================================================================
// vectors for effective rates
//========================================================================
vector<double> Ai_df_dx_HeI, Bi_df_dx_HeI;
vector<vector<double> > RijVec_df_dx_HeI;


//========================================================================
//
// parameters which are reset with parameter file
//
//========================================================================
int nShellsHeI=10;
int njresolved=10;

int flag_HI_absorption=1;         // HI absorption on nP-1s helium lines
                                  // 0: no effect; 
                                  // 1: effect on 2P-1s; 
                                  // 2: including of diffusion correction

int _HI_abs_appr_flag=0;

//========================================================================
// HeI specific switches
//========================================================================
int flag_spin_forbidden=1;
int flag_He=1;

//========================================================================
// print feedback message
//========================================================================
bool He_feedback_message=0;

//========================================================================
// HeI S->T feedback
//========================================================================
int HeISTfeedback=0;
int _HeI_feedback_w_HI_abs=0;
int nHeFeedback=2;                                

//========================================================================
// reference points for rescaling and switch of DPesc
//========================================================================
const double zcrit_HI=3400.0;
const double zcrit_HI_Int=3400.0;
const double zref_HeI=2500.0;

//========================================================================
// both defined in ./Modules/P_HI_abs_tabulation.cpp
//========================================================================
void create_Pesc_tables(double z_s, double z_e, Photoionization_Lyc &Lyc1s, 
                        Gas_of_HeI_Atoms &HeIA, Cosmos &cos);
void read_DP_Data(string filename, int show); 

//========================================================================
// for helium atom
//========================================================================
void setup_Helium_atom()
{
    print_message(" entering the setup for neutral helium ");
    
    if(nShellsHeI<njresolved) njresolved=nShellsHeI;
    
    Set_flag_Intercombination_lines(10);    
    if(set_HeI_Quadrupole_lines)
    {
        Set_flag_Quadrupole_lines(10);
        print_message(" switching on HeI nD-1S singlet quadrupole transitions ");
    }
    
    HeIAtoms=new Gas_of_HeI_Atoms(nShellsHeI, njresolved, -2);
    
    //===========================================================================
    // Load effective rates
    //===========================================================================
    if(nShellsHeI==2) effective_rate_path_HeI=effective_rate_path_HeI_arr[0];
    else if(nShellsHeI==3) 
    {
        if(set_HeI_Quadrupole_lines) effective_rate_path_HeI=effective_rate_path_HeI_arr[4];
        else effective_rate_path_HeI=effective_rate_path_HeI_arr[1];
    }
    else if(nShellsHeI==5)
    {
        if(set_HeI_Quadrupole_lines) effective_rate_path_HeI=effective_rate_path_HeI_arr[5];
        else effective_rate_path_HeI=effective_rate_path_HeI_arr[2];
    }
    else if(nShellsHeI==10)
    {
        if(set_HeI_Quadrupole_lines) effective_rate_path_HeI=effective_rate_path_HeI_arr[6];
        else effective_rate_path_HeI=effective_rate_path_HeI_arr[3];
    }
    else
    { 
        cerr << " setup_Helium_atom :: please choose nSHeI = 2, 3, 5 or 10! Exiting... \n" << endl; 
        exit(0); 
    }
    
    load_rates(effective_rate_path_HeI, nS_effective_HeI, HeI_Atoms);
    
    //===========================================================================
    // allocate memory to pass on the solution to the PDE solver            // JF
    //===========================================================================
    pass_on_the_Solution_CosmoRec_HeI.clear();
    vector<double> dum(1+1+get_number_of_resolved_levels_HeI());
    
    for(int k=0; k<parameters.nz; k++) 
        pass_on_the_Solution_CosmoRec_HeI.push_back(dum);    
    
    //===========================================================================
    // memory for effective rates
    //===========================================================================
    Ai_df_dx_HeI.resize(get_number_of_resolved_levels_HeI());
    Bi_df_dx_HeI.resize(get_number_of_resolved_levels_HeI());
    RijVec_df_dx_HeI.clear();
    for(int k=0; k<(int)Ai_df_dx_HeI.size(); k++) RijVec_df_dx_HeI.push_back(Ai_df_dx_HeI);

    //===========================================================================
    // create/load DP escape tables
    //===========================================================================
    if(flag_HI_absorption!=0) 
        read_DP_Data(Rec_database_path+"/Pesc_Data/DP_Coll_Data.31.fac_50.neff_"
                                      +int_to_string(nS_effective_HeI)+".dat", 1);
    
    return;
}

//========================================================================
void Set_HeI_Levels_to_Saha(double z)
{
    double Xe=cosmos.Xe_Seager(z);
    double Xp=cosmos.Xp(z);
    double NH=cosmos.NH(z);
    double TM=cosmos.TCMB(z);
    double XHeII=min(Xe-Xp, parameters.fHe);
        
    if(show_CosmoRec_mess>=1) 
        cout << " Helium \n" << " Xe-Xp: " << Xe-Xp << " " << XHeII << endl;
    // follow all sub-level 
    for(long int i=0; i<Level_I.nHeIeq; i++)
    {
        Level_I.X[i+Level_I.index_HeI]=HeI_Atoms.Xi_Saha(i, min(Xe, 1.0+XHeII), XHeII, NH, TM); 
        HeI_Atoms.Set_Xi(i, Level_I.X[i+Level_I.index_HeI]);
        
        if(show_CosmoRec_mess>=2) 
        {
            if(i==HeI_Atoms.Get_indexT()) cout << "\n Triplet-levels (j-resolved)" << endl;
            
            if(i==HeI_Atoms.Get_indexT_no_j() && njresolved<nShellsHeI) 
                cout << " Triplet-levels (non-j-resolved)" << endl;
            
            cout << " He-Level: " << i << " " << Level_I.X[i+Level_I.index_HeI] 
                 << " HeI_Atoms: " << HeI_Atoms.Xi(i) << endl;
        }
    }

    if(show_CosmoRec_mess>=1) cout << endl;
    
    return;
}


