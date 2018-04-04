//===========================================================================================================
void print_var_change_message(string mess, double zs)
{
    if(show_CosmoRec_mess<2) return;
    
    cout << "\n " << setfill('*') << setw(100) << "*" << endl;  
    cout << " * " << endl;
    cout << " * " << mess << " " << zs << endl;
    cout << " * " << endl;
    cout << " " << setfill('*') << setw(100) << "*" << endl << endl;   

    return;
}

//===========================================================================================================
void set_accuracies_own_solver(double *rtol, double *abs, int neq_HI, int neq)
{
    for(int k=0; k<neq; k++) rtol[k]=abs[k]=0.0;
    
    rtol[0]=1.0e-7; abs[0]=1.0e-12;
    rtol[1]=1.0e-7; abs[1]=1.0e-12;
    rtol[2]=1.0e-6; abs[2]=1.0e-80;
    rtol[3]=1.0e-6; abs[3]=1.0e-80;
    //
    if(flag_He==1)
    {
        rtol[neq_HI+1]=1.0e-7; abs[neq_HI+1]=1.0e-12;
        rtol[neq_HI+2]=1.0e-6; abs[neq_HI+2]=1.0e-80;
        rtol[neq_HI+3]=1.0e-6; abs[neq_HI+3]=1.0e-80;
    }
    
    return;
}

//===========================================================================================================
//
// The calculation
//
//===========================================================================================================
int Xe_frac_effective_rates()
{
    bool numerical_jacobian=0;
    int err=0;

    //==========================================================================
    // HeI S->T feedback tables
    //==========================================================================
    if(HeISTfeedback==1) init_HeI_feedback(HeI_Atoms, parameters, zarr);
    
    //===========================================================================
    // initial values for populations (Level I)
    //===========================================================================
    if(show_CosmoRec_mess>=1) print_message(" Setting initial populations... "); 
    
    //===========================================================================
    // hydrogen & helium specific setup
    //===========================================================================
    Set_Hydrogen_Levels_to_Saha(parameters.zstart);
    Set_HeI_Levels_to_Saha(parameters.zstart);
    
    //===========================================================================
    // Now He I recombination starts
    //===========================================================================
    if (parameters.zstart > parameters.zend) 
    {
        double zs=zarr[0], ze=0.0;
        
        //===========================================================================
        // variables for solver
        //===========================================================================
        int neq=3+get_number_of_resolved_levels()+get_number_of_resolved_levels_HeI(); 
        
        //===========================================================================
        // set up Solution memory for my ODE-solver
        //===========================================================================
        ODE_solver_Rec::allocate_memory(Sz, tols, neq);
        Sz.z=zs;        
        set_accuracies_own_solver(&tols.rel[0], &tols.abs[0], get_number_of_resolved_levels(), neq);
        
        //===========================================================================
        // initialize Sz
        //===========================================================================
        initialize_y(Level_I, Hydrogen_Atoms, HeI_Atoms, zs, Level_III);
        copy_LIII_to_ysol(Level_III.y, &Sz.y[0], Level_III.neq, Level_III.index_HI, Level_III.index_HeI);
        
        //===========================================================================
        // initialize ODE solver
        //===========================================================================
        ODE_solver_Rec::ODE_Solver_set_up_solution_and_memory(Sz, tols, ODE_Solver_info_CR, 
                                                              fcn_effective, numerical_jacobian);
        
        //===========================================================================
        // files etc
        //===========================================================================
        open_files(Xfiles, path);   

        //===========================================================================
        // general output
        //===========================================================================
        write_files(zs, Xfiles);
        save_sol_for_PDE_solver(0, zs);
        save_Xe_Te_for_output(zs);
        
        if(iteration_count==0) print_message(" Entering main ODE computation. ");
        else  print_message(" Re-entering main ODE computation. ");
        
        for(int i=1; i<parameters.nz; i++)
        {
            zs=zarr[i-1]; ze=zarr[i];
            
            if( (zs<200.0 || parameters.fHe-HeI_Atoms.Xi(0)<=Xi_HeI_switch) && flag_He==1)
            { 
                print_var_change_message("Xe_frac :: switching off helium at XHeII =", parameters.fHe-HeI_Atoms.Xi(0));
                flag_He=0;
                HeI_was_switched_off_at_z=zs;
                
                //--------------------------------------------------------
                // set helium to be completely recombined
                //--------------------------------------------------------
                Level_I.X[Level_I.index_HeI]=parameters.fHe;
                for(int i=Level_I.index_HeI+1; i<Level_I.index_HeI+Level_I.nHeIeq; i++) Level_I.X[i]=0.0;
                for(int i=0; i<(int)HeI_Atoms.Get_total_number_of_Levels(); i++) HeI_Atoms.Set_Xi(i, Level_I.X[Level_I.index_HeI+i]);
                
                //============================================
                // take away helium equations
                //============================================
                neq=2+get_number_of_resolved_levels();
                set_variables(700.0, zref_HeI);
                initialize_y(Level_I, Hydrogen_Atoms, HeI_Atoms, zs, Level_III);
                
                //============================================
                ODE_solver_Rec::allocate_memory(Sz, tols, neq);
                Sz.z=zs;        
                set_accuracies_own_solver(&tols.rel[0], &tols.abs[0], get_number_of_resolved_levels(), neq);
                
                copy_LIII_to_ysol(Level_III.y, &Sz.y[0], Level_III.neq, Level_III.index_HI, Level_III.index_HeI);
                
                ODE_solver_Rec::ODE_Solver_set_up_solution_and_memory(Sz, tols, ODE_Solver_info_CR, fcn_effective, numerical_jacobian);
            }
            
            //===========================================================================
            // main run
            //===========================================================================
            err=ODE_solver_Rec::ODE_Solver_Solve_history(zs, ze, Sz, ODE_Solver_info_CR);            
            if(err==10) break;            
            copy_ysol_to_LIII(&Sz.y[0], Level_III.y, Level_III.neq, Level_III.index_HI, Level_III.index_HeI);
            
            //===========================================================================
            // save step
            //===========================================================================
            copy_y_to_X(Level_III, ze, Level_I);
            //
            Level_I.X[0]=Level_III.y[0]=1.0-Level_I.X[Level_I.index_HI];
            update_H_populations(&Level_I.X[Level_I.index_HI], Hydrogen_Atoms);
            if(flag_He==1)
            {
                update_HeI_populations(&Level_I.X[Level_I.index_HeI], HeI_Atoms);
                Level_I.X[0]+= parameters.fHe-Level_I.X[Level_I.index_HeI];
                Level_III.y[0]=Level_I.X[0];
            }
            
            //===========================================================================
            // outputs
            //===========================================================================
            write_files(ze, Xfiles);
            
            //===========================================================================
            // Save results for PDE solver
            //===========================================================================
            save_sol_for_PDE_solver(i, ze);
            
            //===========================================================================
            // Save results for output to CosmoMC
            //===========================================================================
            save_Xe_Te_for_output(ze);
            
            if(show_CosmoRec_mess>=1)
            {
                cout << endl;
                cout << "\n *----------------------------------------------------------------------------------------------------------* " << endl;
                cout << " current redshift: i= " << i << " z= " << ze << " Xe= " << Level_I.X[0] << " XHI1s= " << Level_I.X[Level_I.index_HI];
                if(flag_He==1) cout << " XHeI1s= " << Level_I.X[Level_I.index_HeI];
                cout << " rho= " << Level_I.X[Level_I.neq-1] << endl;
            }
            
            //==========================================================================
            // HeI S->T feedback (some messages can be printed here...)
            //==========================================================================
            // 22.02.2011: this function is not evalutated when the helium diffusion
            // correction is treated with full PDE-solver.
            //==========================================================================
            if(HeISTfeedback==1 && flag_He==1 && !Diffusion_correction_HeI_is_on) 
                save_HeI_distortion(HeI_Atoms, Hydrogen_Atoms.X(0), i, ze);
            
            if(show_CosmoRec_mess>=1)
                cout << " *----------------------------------------------------------------------------------------------------------* " << endl << endl;
        }
        
        if(flag_compute_Recfast_part==1 && err==0) 
        {
            //==========================================================================
            // set derivative of Xe for rescaling 
            //==========================================================================
            fcn_effective(ze, Level_III);
            double dXe_dz=-Level_III.g[Level_I.index_HI]*f0[Level_I.index_HI]/(-cosmos.H(ze)*(1.0+(ze))); 
            
            compute_Recfast_part(ze, ze_Recfast, 1.0-Hydrogen_Atoms.X(0), 0.0, Level_I.X[0], 
                                 Level_I.X[Level_I.neq-1], dXe_dz, Xfiles[1], Xfiles[0]);
        }
        
        //==============================================================================
        // Clean up  
        //==============================================================================
        close_files(Xfiles);            
    }
    
    return err;
}

//======================================================================================
//
// Main call
//
//======================================================================================
int compute_history_with_effective_rates()
{ return Xe_frac_effective_rates(); }

