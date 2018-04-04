//===========================================================================================================
// Author Jens Chluba 
// here functions which are common for all the programs
// connected with the recombination problem are defined
//===========================================================================================================

//===========================================================================================================
void set_startup_data(double params[12], Parameters_form_initfile &inparams)
{
    //=======================================================================================================
    // params has to contain 11 parameters in the following sequence: 
    // nzpts, zs, ze, Y_p, T0, OmT, OmB, OmL, OmK, h100, fudge
    // if the fudge factor ==0.0 --> will use standard value 1.14
    //=======================================================================================================
    
    inparams.nz=(int)params[0];
    inparams.zstart=min(3400.0, params[1]);
    inparams.zend=params[2];
    inparams.YP=params[3];
    inparams.To=params[4];
    inparams.Omega0=params[5];
    inparams.OmegaB=params[6];
    inparams.OmegaL=params[7];
    inparams.OmegaK=params[8];
    inparams.h100 = params[9];
    inparams.Nnu = params[10];
    if(params[11]==0.0) inparams.F=1.14;
    else    inparams.F=params[11];
    
    inparams.H0 = params[9]*100.0*1.0e+5/const_Mpc;
    inparams.fHe=cosmos.fHe(params[3]);
    
    return;
}

//===========================================================================================================
void set_cosmology(Parameters_form_initfile &inparams, Cosmos &cosmos, double fDM_ann, int CR_corr)
{
    //=======================================================================================================
    // for the interpolation function for Xe in Cosmos class
    //=======================================================================================================
    int n_Xe=10000;
    double dens[7]; 
    //=======================================================================================================
    
    if(show_CosmoRec_mess>=1) print_message(" Initializing cosmology object ");
    
    dens[0]=inparams.Omega0+inparams.OmegaL+inparams.OmegaK; 
    dens[1]=inparams.Omega0*pow(inparams.h100, 2);
    dens[2]=inparams.OmegaB*pow(inparams.h100, 2);
    dens[3]=inparams.OmegaL;
    dens[4]=inparams.F;
    dens[5]=fDM_ann;
    dens[6]=CR_corr;
    
    //=======================================================================================================
    // initializing cosmos object...
    //=======================================================================================================
    cosmos.init(inparams.h100, inparams.To, inparams.YP, dens, inparams.zstart, 1.0e-10, n_Xe, inparams.Nnu);
    cosmos.init_splines();
    
    if(show_CosmoRec_mess>=1) cosmos.display_cosmology();
    
    return;
}

//===========================================================================================================
// initialization
//===========================================================================================================
void set_array_dimensions_and_parameters(int wait=0)
{
    //=======================================================================================================
    // setup Level I
    //=======================================================================================================
    Level_I.nHIeq=Hydrogen_Atoms.Get_total_number_of_Levels();
    Level_I_temp.nHIeq=Level_I.nHIeq;
    //
    Level_I.nHeIeq=HeI_Atoms.Get_total_number_of_Levels();
    Level_I_temp.nHeIeq=Level_I.nHeIeq;
    //
    Level_I.neq=Level_I.nHIeq+Level_I.nHeIeq+2;
    Level_I_temp.neq=Level_I.neq;
    
    //=======================================================================================================
    // start index of HI and HeI
    //=======================================================================================================
    Level_I.index_HI=1;
    Level_I_temp.index_HI=Level_I.index_HI;
    //
    Level_I.index_HeI=Level_I.nHIeq+1;
    Level_I_temp.index_HeI=Level_I.index_HeI;
    
    //=======================================================================================================
    // setup Level III
    //=======================================================================================================
    Level_III.nHIeq=Level_III_temp.nHIeq=Level_III_temp_JAC.nHIeq=Hydrogen_Atoms.Get_number_of_Levels_until(nShells);
    Level_III.nHeIeq=Level_III_temp.nHeIeq=Level_III_temp_JAC.nHeIeq=HeI_Atoms.Get_total_number_of_Levels();
    Level_III.neq=Level_III_temp.neq=Level_III_temp_JAC.neq=Level_III.nHIeq+Level_I.nHeIeq+2;

    //=======================================================================================================
    // start index of HI and HeI (these are reset when using reordering)
    //=======================================================================================================
    Level_III.index_HI =Level_III_temp.index_HI=Level_III_temp_JAC.index_HI=1;
    Level_III.index_HeI=Level_III_temp.index_HeI=Level_III_temp_JAC.index_HeI=Level_III.nHIeq+1;

    //=======================================================================================================
    if(show_CosmoRec_mess>=2)
    {
        cout << "\n total number of equations: " << Level_III.neq 
             << "\n Hydrogen equations: " << Level_III.nHIeq       
             << "\n Helium equations: " << Level_III.nHeIeq << endl;
    
        cout << " nShells: " << nShells << "\n included number of hydrogen levels: " 
             << Level_I.nHIeq << endl;
    
        cout << " nShellsHeI: " << nShellsHeI 
             << "\n included number of neutral helium levels: " 
             << Level_I.nHeIeq << endl;
    }
    //=======================================================================================================

    if(wait==1) wait_f_r();
    
    return;
}

//===========================================================================================================
void read_entries_from_parameter_file(ifstream &paramfile, double *param)
{
    for(int i=0; i<12; i++) paramfile >> param[i];

    paramfile >> nShells;
    paramfile >> nS_effective;
    paramfile >> fDM_CosmoRec;
    
    //======================================================================
    // explicitly set value of A2s1s rate (added 23.07.2014)
    //======================================================================
    double dumA2s1s;
    paramfile >> dumA2s1s;
    if(dumA2s1s>0.0) const_HI_A2s_1s=RECFAST_atomic_constants::RF_Lam2s1sH=dumA2s1s;

    paramfile >> nShellsHeI;
    paramfile >> flag_HI_absorption;
    paramfile >> flag_spin_forbidden;
    paramfile >> nHeFeedback;
    
    paramfile >> Diffusion_flag;
    paramfile >> induced_flag;
    paramfile >> ntg_max;
    paramfile >> nR_max;    
        
    paramfile >> path;
    paramfile >> addname;

    return;
}

//===========================================================================================================
void set_detailed_name();                                                    // defined in output_methods.cpp

//===========================================================================================================
//
// setting the startup parameters according to param[12]
//
//===========================================================================================================
void set_startup_data(string filename, double param[12], Parameters_form_initfile &inparams, 
                      Cosmos &cos, int wait=0)
{
    if(show_CosmoRec_mess>=0) cout << endl;
    
    //=======================================================================================================
    // compute Omega_L according to Omega_K
    //=======================================================================================================
    if(param[7]<=0.0) param[7]=1.0-param[5]-param[8]-cos.calc_Orel(param[4], param[10], param[9]);
    
    if(nHeFeedback==0 || nHeFeedback==1 || nHeFeedback==-1)
    { 
        nHeFeedback=0; 
        HeISTfeedback=0; 
        _HeI_feedback_w_HI_abs=0; 
    }
    else if(nHeFeedback>0)
    { 
        nHeFeedback=(int)min(nHeFeedback, nShellsHeI); 
        HeISTfeedback=1; 
        _HeI_feedback_w_HI_abs=0; 
    }
    else 
    { 
        nHeFeedback=(int)min(-nHeFeedback, nShellsHeI); 
        HeISTfeedback=1; 
        _HeI_feedback_w_HI_abs=1; 
    }
    
    if(flag_HI_absorption==1){}
    else if(flag_HI_absorption==2){ _HI_abs_appr_flag=1; flag_HI_absorption=1; }
    else if(flag_HI_absorption==3){ _HI_abs_appr_flag=1; flag_HI_absorption=1; Diffusion_correction_HeI=1; }
    else if(flag_HI_absorption>=4)
    { 
        cerr << " set_startup_data :: flag_HI_absorption flag N/A. Exiting... " << endl; 
        exit(0);
    }

    if(fDM_CosmoRec>0) DM_annihilation=1;
    
    if(show_CosmoRec_mess>=1)
    {
        if(filename!="") cout << " Using parameters corresponding to " << filename << endl;
        cout << "\n zs: " << param[1] << "\t ze: " << param[2] << "\t nzpts: " 
             << (int)param[0] << endl;
        cout << " Y_p: " << param[3] << "\t TCMB: " << param[4] << endl;
        cout << " OmegaM: " << param[5] << "\t OmegaB: " << param[6] 
             << "\n OmegaL: " << param[7] << "\t OmegaK: " << param[8] 
             << "\t Omega_rel: " << cos.calc_Orel(param[4], param[10], param[9]) 
             << " Nnu= " << param[10] << endl;
        cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << param[9] << endl;
        cout << " Fudge factor for H recombination: " << param[11] << endl;
        cout << " HI A2s1s rate: " << const_HI_A2s_1s << endl;
    }
    
    if(DM_annihilation && show_CosmoRec_mess>=0)
        cout << " Including dark matter annihilation. f_DM= " << fDM_CosmoRec << " eV/sec " << endl;
    
    if(flag_HI_absorption==1 && show_CosmoRec_mess>=0) 
            cout << " HI absorption will be included during HeI-recombination " << endl;
    //
    if(Diffusion_correction_HeI==1 && show_CosmoRec_mess>=0) 
        cout << " After the first iteration of the normal EMLA, the HeI-diffusion" 
             << " correction will be computed from PDE-solver " << endl;
    //
    if(flag_spin_forbidden>=1)
    { 
        if(show_CosmoRec_mess>=0) 
            cout << " HeI 2P-1S intercombination line is switched on " << endl; 
        flag_spin_forbidden=1; 
    }
    else
    { 
        if(show_CosmoRec_mess>=0) 
            cout << " HeI 2P-1S intercombination line is switched off " << endl; 
        flag_spin_forbidden=0; 
    }

    if(Atom_activate_HI_Quadrupole_lines)
    { 
        if(show_CosmoRec_mess>=0 && nShells==3) 
            cout << " HI quadrupole lines are included " << endl; 
    }
    
    if(Diffusion_flag>0)
    { 
        Diffusion_flag=1;
        Diffusion_correction=1;
        DI1_2s_correction_on=0;
        switch_off_1s_2s_correction();
        
        if(show_CosmoRec_mess>=0)
        {
            cout << " Diffusion correction for HI Ly-a will be included." 
                 << " The number of iterations will be "+int_to_string(Diff_iteration_max) 
                 << endl;
        }

        if(induced_flag>=1)
        { 
            switch_on_1s_2s_correction();
            DI1_2s_correction_on=1;
            
            if(show_CosmoRec_mess>=0) 
                cout << " Switching on stimulated 2s-1s two-photon decays " << endl; 
            
            if(induced_flag==2)
            { 
                if(show_CosmoRec_mess>=0)
                    cout << " Ly-n feedback to the 1s-2s two-photon channel " << endl; 
            }
            else switch_off_1s_2s_absorption_correction();     // defined in Solve_PDEs.h
        }
    }
        
    //=======================================================================================================
    // below given z extend the solution not solving the whole system
    //=======================================================================================================
    if(param[2]<50.0)
    {
        flag_compute_Recfast_part=1;
        nz_Recfast_part=100;
        ze_Recfast=max(param[2], 0.001);
        param[2]=50.0;
    }   
    
    set_startup_data(param, inparams);
    
    if( show_CosmoRec_mess>=0 && (write_Xe_Te_sol!=0 || write_populations!=0) )
    {
        cout << "\n output-path: " << path << endl;
        cout << " addname    : " << addname << endl;
    }
    
    set_detailed_name();
    
    if(show_CosmoRec_mess>=0) cout << endl;
    
    if(wait==1) wait_f_r();
    
    set_cosmology(inparams, cos, 0.0, 0);       

    if(!hydrogen_and_helium_atom_are_set)
    { 
        setup_Hydrogen_atom();
        setup_Helium_atom();
        hydrogen_and_helium_atom_are_set=1;
    }
    
    return;
}

//===========================================================================================================
// reading startup parameters from startup file
//===========================================================================================================
void read_startup_data(string filename, Parameters_form_initfile &inparams, 
                       Cosmos &cos, int wait=0)
{
    ifstream paramfile(filename.c_str());
    if(!paramfile)
    { 
        cerr << " read_startup_data :: Error opening parameter file. Exiting. " << endl; 
        exit(0); 
    }
    
    double param[12];
    read_entries_from_parameter_file(paramfile, param);
    paramfile.close();
    
    set_startup_data(filename, param, inparams, cos, wait);
    
    return;
}   

//===========================================================================================================
// reset to fresh startup
//===========================================================================================================
void restart_CosmoRec()
{   
    flag_He=1;
    set_array_dimensions_and_parameters();
    set_variables();
    
    return;
}   

//========================================================================================
// Save results for PDE solver
//========================================================================================
void save_sol_for_PDE_solver(int iz, double z)
{
    if(!Diffusion_correction || (Diffusion_correction && iteration_count<Diff_iteration_max) )
    {
        pass_on_the_Solution_CosmoRec[iz][0]=z;
        
        for(int l=0; l<Level_I.nHIeq+1; l++) 
            pass_on_the_Solution_CosmoRec[iz][l+1]=Level_I.X[l];

        pass_on_the_Solution_CosmoRec[iz][Level_I.nHIeq+2]=Level_I.X[Level_I.neq-1];
        pass_on_the_Solution_CosmoRec[iz][Level_I.nHIeq+3]=0.0;
        pass_on_the_Solution_CosmoRec[iz][Level_I.nHIeq+4]=0.0;   

        //================================================================================
        // Helium diffusion correction; only resolved states + XHeI1s are saved      // JF
        //================================================================================
        if(Diffusion_correction_HeI)
        {
            pass_on_the_Solution_CosmoRec_HeI[iz][0]=z;
            pass_on_the_Solution_CosmoRec_HeI[iz][1]=Level_I.X[Level_I.index_HeI];
            
            for(int l=0; l<get_number_of_resolved_levels_HeI(); l++) 
                pass_on_the_Solution_CosmoRec_HeI[iz][l+2]=Level_I.X[Level_I.index_HeI+get_HeI_index(l)];
        }
    }
    
    return;
}

//========================================================================================
// Save results for output to CosmoMC
//========================================================================================
void save_Xe_Te_for_output(double z)
{
    if(!Diffusion_correction || (Diffusion_correction && iteration_count==Diff_iteration_max) )
    {
        vector<double> dum(3);
        dum[0]=z; 
        dum[1]=Level_I.X[0]; 
        dum[2]=Level_I.X[Level_I.neq-1]*cosmos.TCMB(z);
        
        output_CosmoRec.push_back(dum);
    }

    return;
}

//===========================================================================================================
// memory (de)allocation
//===========================================================================================================
void clear_atoms()
{
    delete HAtoms;
    delete HeIAtoms;
    
    Ai_df_dx.clear();
    Bi_df_dx.clear();
    RijVec_df_dx.clear();
    //
    Ai_df_dx_HeI.clear();
    Bi_df_dx_HeI.clear();
    RijVec_df_dx_HeI.clear();
    
    output_CosmoRec.clear();
    pass_on_the_Solution_CosmoRec.clear();
    pass_on_the_Solution_CosmoRec_HeI.clear(); // JF
    
    return;
}

//===========================================================================================================
void allocate()
{
    zarr =new double[parameters.nz];
    // Level I
    Level_I.X =new double[Level_I.neq];
    Level_I.g =new double[Level_I.neq];
    //
    Level_I_temp.X =new double[Level_I.neq];
    Level_I_temp.g =new double[Level_I.neq];
    
    // Level III
    Level_III.y   =new double[Level_III.neq];
    Level_III.g   =new double[Level_III.neq];
    //
    Level_III_temp.y   =new double[Level_III.neq];
    Level_III_temp.g   =new double[Level_III.neq];
    //
    Level_III_temp_JAC.y   =new double[Level_III.neq];
    Level_III_temp_JAC.g   =new double[Level_III.neq];
    //----------------------------------------------------
    
    f0=new double[Level_I.neq];
    
    ODE_Solver_info_CR.Snewptr=NULL; // to make sure that the memory will be set

    return;
}

//===========================================================================================================
void deallocate()
{
    delete[] zarr;
    // Level I
    delete[] Level_I.X;
    delete[] Level_I.g;
    //
    delete[] Level_I_temp.X;
    delete[] Level_I_temp.g;
    
    // Level III
    delete[] Level_III.y;
    delete[] Level_III.g;
    //
    delete[] Level_III_temp.y;
    delete[] Level_III_temp.g;
    //
    delete[] Level_III_temp_JAC.y;
    delete[] Level_III_temp_JAC.g;
    
    delete[] f0;
    delete[] Xfiles;
    
    clear_all_memory(Sz, tols, ODE_Solver_info_CR);
    
    return;
}

//===========================================================================================================
//===========================================================================================================
