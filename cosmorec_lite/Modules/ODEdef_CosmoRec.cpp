//==================================================================================================
// Hydrogen atom
//==================================================================================================
void update_H_populations(const double *XiH, Gas_of_Atoms &H_Atoms)
{
    // XiH[0..nH] contains ALL the populations of levels
    H_Atoms.Set_population_of_level(0, XiH[0]);
    for(int i=0; i<get_number_of_resolved_levels(); i++) 
        H_Atoms.Set_population_of_level(get_HI_index(i), XiH[get_HI_index(i)]);
    
    return;
}

//==================================================================================================
// neutral Helium setup
//==================================================================================================
void update_HeI_populations(const double *XiHeI, Gas_of_HeI_Atoms &HeIA)
{
    // XiHeI[0..nH] contains ALL the populations of levels
    HeIA.Set_Xi(0, XiHeI[0]);
    for(int i=0; i<get_number_of_resolved_levels_HeI(); i++) 
        HeIA.Set_Xi(get_HeI_index(i), XiHeI[get_HeI_index(i)]);
    
    return;
}

//==================================================================================================
//==================================================================================================
static double ODE_Solver_z_pre=0.0, ODE_Solver_rho_pre=0.0;

//==================================================================================================
void update_Solver_vars_eff(const double *XiH, Gas_of_Atoms &H_Atoms, 
                            const double *XiHeI, Gas_of_HeI_Atoms &HeIA, 
                            double z, double rho, int Pop_flag=1)
{
    // update populations
    if(Pop_flag==1)
    {
        update_H_populations(XiH, H_Atoms);
        if(flag_He==1) update_HeI_populations(XiHeI, HeIA);
    }
    
    if(z!=ODE_Solver_z_pre || rho!=ODE_Solver_rho_pre)
    {
        ODE_Solver_z_pre=z;
        ODE_Solver_rho_pre=rho;
        
        // effective rates
        get_rates(cosmos.TCMB(z), cosmos.TCMB(z)*rho, Ai_df_dx, Bi_df_dx, RijVec_df_dx);
        if(flag_He==1) get_rates_HeI(cosmos.TCMB(z), Ai_df_dx_HeI, Bi_df_dx_HeI, RijVec_df_dx_HeI);
    }
    
    return;
}

//==================================================================================================
// aux-functions for variable mapping
//==================================================================================================
void copy_LIII_to_ysol(double *LIIIy, double *ysol, int neq, int index_HI, int index_HeI)
{
    int nres_HI=get_number_of_resolved_levels();
    int nres_HeI=get_number_of_resolved_levels_HeI();
    
    //============================================================================
    // 1. eq is Te
    //============================================================================
    ysol[0]=LIIIy[neq-1];
    ysol[1]=LIIIy[index_HI];
    //
    for(int k=0; k<nres_HI; k++){ ysol[k+2]=LIIIy[index_HI+get_HI_index(k)]; }
    //
    if(flag_He==1) 
    {
        ysol[nres_HI+2]=LIIIy[index_HeI]; 
        for(int k=0; k<nres_HeI; k++){ ysol[nres_HI+k+3]=LIIIy[index_HeI+get_HeI_index(k)]; }
    }
    
    return;
}

//==================================================================================================
void copy_ysol_to_LIII(double *ysol, double *LIIIy, int neq, int index_HI, int index_HeI)
{
    int nres_HI=get_number_of_resolved_levels();
    int nres_HeI=get_number_of_resolved_levels_HeI();
    
    //============================================================================
    // 1. eq is Te
    //============================================================================
    LIIIy[neq-1]=ysol[0];
    LIIIy[index_HI]=ysol[1];
    //  
    for(int k=0; k<nres_HI; k++){ LIIIy[index_HI+get_HI_index(k)]=ysol[k+2]; }
    //
    if(flag_He==1) 
    {
        LIIIy[index_HeI]=ysol[nres_HI+2];
        for(int k=0; k<nres_HeI; k++){ LIIIy[index_HeI+get_HeI_index(k)]=ysol[nres_HI+k+3]; }
    }
    
    return;
}

//==================================================================================================
// fcn for helium
//==================================================================================================
void fcn_He_effective(double z, double Xe, double Tg, double NH, double H_z, 
                      double XHeII, int iHeI, Data_Level_I &LI)
{   
    double Te=Tg;   
    
    //==================================================================
    // these are the ones with effective rates that depend on Xi
    //==================================================================    
    ODE_HeI_effective::evaluate_effective_Rci_Ric_terms(z, Xe, XHeII*NH, HeI_Atoms, Ai_df_dx_HeI, 
                                                        Bi_df_dx_HeI, get_HeI_index, &LI.g[iHeI]);
    
    ODE_HeI_effective::evaluate_effective_Rij_terms(z, Tg, HeI_Atoms, RijVec_df_dx_HeI, 
                                                    get_HeI_index, &LI.g[iHeI]);    
    
    //==========================================================================
    // 2s two-photon channels (2S - 1S singlet channel)
    //==========================================================================
    ODE_effective::evaluate_2s_two_photon_decay(z, Tg, HeI_Atoms.Xi(0), HeI_Atoms.Xi(1), 
                                                HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                                const_HeI_A2s_1s, LI.g[iHeI], LI.g[iHeI+1]);    
    
    //==================================================================
    // HeI 'Lyman-series' & quadrupole lines (if switched on)
    //==================================================================
    for(int k=1; k<get_number_of_resolved_levels_HeI(); k++)
    {
        int ik=get_HeI_index(k);
        if(flag_spin_forbidden==0 && ik==(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1)) continue;
        
        Transition_Data_HeI_A T=HeI_Atoms.Get_Trans_Data(ik, 1, 0, 0, 0);
        if(T.A21!=0) ODE_effective::evaluate_Ly_n_channel(z, Tg, HeI_Atoms.Xi(0), HeI_Atoms.Xi(ik), 
                                                          NH, H_z, T.A21, T.lambda21, T.Dnu, 
                                                          LI.g[iHeI], LI.g[iHeI+ik], 
                                                          HeI_Atoms.Get_gw(ik)/T.gwp);
    }
    
    //==================================================================
    // HeI-diffusion correction (added 22.01.2011; JF & JC)
    //==================================================================
    if(Diffusion_correction_HeI_is_on)
    {
        evaluate_HeI_Diffusion_correction(z, HeI_Atoms, NH, H_z, Tg, 
                                          &LI.g[iHeI], LI.g[0], LI.g[LI.index_HI]);
        
        if(DI1_2s_correction_on_HeI) 
            evaluate_HeI_DI1_DI2_correction(z, HeI_Atoms.Xi(0), HeI_Atoms.Xi(1), 
                                            HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                            Tg, &LI.g[iHeI]);
    }
    
    //==================================================================
    // HeI-feedback
    //==================================================================
    if(HeISTfeedback==1 && !Diffusion_correction_HeI_is_on)
    {
        if(_HeI_feedback_w_HI_abs==0) 
            evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                  Dn_HeI_feedback_no_HI_abs, &LI.g[iHeI]);
        
        if(_HeI_feedback_w_HI_abs==1) 
            evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                  Dn_HeI_feedback, &LI.g[iHeI]);
    }
    
    //==================================================================
    // absorption by neutral hydrogen
    //==================================================================
    if(flag_HI_absorption==1 && z<=zcrit_HI && !Diffusion_correction_HeI_is_on)
    {  
        double dXe_HeI_HI_abs_dt=0.0;
        //==============================================================
        // Ly-a line
        //==============================================================
        evaluate_HI_abs_HeI(z, Hydrogen_Atoms.X(1, 0), NH, H_z, Tg, Te, 
                            HeI_Atoms, cosmos, 
                            HeI_Atoms.HILyc, DP_interpol_S, dXe_HeI_HI_abs_dt);

        // Helium levels
        LI.g[iHeI+0]+= dXe_HeI_HI_abs_dt;
        LI.g[iHeI+2]+=-dXe_HeI_HI_abs_dt;
        // electron equation
        LI.g[0]+=dXe_HeI_HI_abs_dt;
        // HI 1s
        LI.g[LI.index_HI]+=-dXe_HeI_HI_abs_dt;
        
        //==============================================================
        // from intercombination line
        //==============================================================
        if(z<=zcrit_HI_Int && flag_spin_forbidden==1)
        {
            dXe_HeI_HI_abs_dt=0.0;
            evaluate_HI_abs_HeI_Intercombination(z, Hydrogen_Atoms.X(1, 0), NH, H_z, Tg, Te,
                                                 HeI_Atoms, cosmos, 
                                                 HeI_Atoms.HILyc, DP_interpol_T, 
                                                 dXe_HeI_HI_abs_dt);        
            
            // Helium levels
            LI.g[iHeI+0]+= dXe_HeI_HI_abs_dt;
            LI.g[iHeI+HeI_Atoms.Get_Level_index(2, 1, 1, 1)]+=-dXe_HeI_HI_abs_dt;
            // electron equation
            LI.g[0]+=dXe_HeI_HI_abs_dt;
            // HI 1s
            LI.g[LI.index_HI]+=-dXe_HeI_HI_abs_dt;
        }
    }
    
    return;
}

//==================================================================================================
// fcn for ODE evaluation
//==================================================================================================
void fcn_effective(double z, Data_Level_III &LIII)
{
    //==================================================================
    // transform from yi --> Xi
    //==================================================================
    transform_y2X(LIII, z, Level_I_temp);
    
    double rho=Level_I_temp.X[Level_I_temp.neq-1];
    double Tg=cosmos.TCMB(z);
    int iHI=Level_I_temp.index_HI;
    int iHeI=Level_I_temp.index_HeI;
    //
    double Xp=1.0-Level_I_temp.X[iHI];
    double XHeII = (flag_He==1 ? parameters.fHe-Level_I_temp.X[iHeI] : 0.0); 
    double Xe=Xp+XHeII;
    double NH=cosmos.NH(z), H_z=cosmos.H(z), Np=NH*Xp;
    
    //==================================================================
    // copy Xi_H --> Atoms && update recombination coefficients
    //==================================================================
    update_Solver_vars_eff(&Level_I_temp.X[iHI], Hydrogen_Atoms, 
                           &Level_I_temp.X[Level_I_temp.index_HeI], 
                           HeI_Atoms, z, rho);
    
    //==================================================================
    // reset everything 
    //==================================================================
    for(int i=0; i<Level_I_temp.neq; i++) Level_I_temp.g[i]=LIII.g[i]=0.0;        
    
    //==================================================================
    // effective rates
    //==================================================================
    ODE_effective::evaluate_TM(z, Xe, parameters.fHe, rho, Tg, H_z, 
                               Level_I_temp.g[Level_I_temp.neq-1]);
    
    ODE_effective::evaluate_2s_two_photon_decay(z, Tg, Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(1), 
                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                const_HI_A2s_1s, Level_I_temp.g[iHI], 
                                                Level_I_temp.g[iHI+1]);
    
    //==================================================================
    // HI Lyman-series
    //==================================================================
    int nres=get_number_of_resolved_levels();
    int nLy_max=(int)min(Lyn_max_parameter, Hydrogen_Atoms.Get_n_of_Level(get_HI_index(nres-1)));
    
    for(int n=2; n<=nLy_max; n++)
    {
        Transition_Data_A T=Hydrogen_Atoms.Level(n, 1).Get_Data_A(1, 0);
        
        ODE_effective::evaluate_Ly_n_channel(z, Tg, Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(n, 1), 
                                             NH, H_z, T.A21, T.lambda21, T.Dnu, 
                                             Level_I_temp.g[iHI], 
                                             Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n,1)]);
    }   
    
    //==================================================================
    // HI nD-1s Quadrupole transitions (added May 2011)
    //==================================================================
    if(Atom_activate_HI_Quadrupole_lines)
    {
        for(int k=0; k<get_number_of_resolved_levels(); k++)
        {
            int l=Hydrogen_Atoms.Get_l_of_Level(get_HI_index(k));
            
            if(l==2)
            {
                int n=Hydrogen_Atoms.Get_n_of_Level(get_HI_index(k));
                
                if(n<=nD_Quadrupole_max)
                    ODE_effective::evaluate_nD_1s_Q_channel(Hydrogen_Atoms.X(n, 1), 
                                                            Hydrogen_Atoms.X(n, 2), 
                                                            Hydrogen_Atoms.Level(n, 2).Get_A21(1, 0), 
                                                            Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 1)], 
                                                            Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 2)]);
            }
        }
    }

    //==================================================================
    // these are the ones with effective rates that depend on Te
    //==================================================================    
    ODE_HI_effective::evaluate_effective_Rci_Ric_terms(z, Xe, Np, Hydrogen_Atoms, Ai_df_dx, 
                                                       Bi_df_dx, get_HI_index, &Level_I_temp.g[iHI]);
    
    ODE_HI_effective::evaluate_effective_Rij_terms(z, Tg, Hydrogen_Atoms, RijVec_df_dx, 
                                                   get_HI_index, &Level_I_temp.g[iHI]);    
    
    //==================================================================
    // Diffusion correction
    //==================================================================
    if(Diffusion_correction_is_on)
    {
        evaluate_HI_Diffusion_correction(z, Hydrogen_Atoms, NH, H_z, Tg, &Level_I_temp.g[iHI]);
        
        if(DI1_2s_correction_on) evaluate_HI_DI1_DI2_correction(z, Hydrogen_Atoms.X(0), 
                                                                Hydrogen_Atoms.X(1), 
                                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                                Tg, &Level_I_temp.g[iHI]);
    }
    
    //==================================================================
    // DM annihilation
    //==================================================================
    if(DM_annihilation) 
        evaluate_DM_annihilation_terms(z, parameters.fHe, Tg, Xp, XHeII, fDM_CosmoRec, 
                                       Level_I_temp.g[iHI], Level_I_temp.g[iHeI], 
                                       Level_I_temp.g[Level_I_temp.neq-1]);
        
    //==================================================================
    // neutral helium part
    //==================================================================
    if(flag_He==1) fcn_He_effective(z, Xe, Tg, NH, H_z, XHeII, iHeI, Level_I_temp);
    
    //==================================================================
    // to transform g_i --> dy_i/dt
    //==================================================================
    transform_g_tX2g_ty(Level_I_temp, z, LIII);
    
    return;
}

//==================================================================================================
//
// jac for ODE evaluation
//
//==================================================================================================
double df_dx_2point_sym(const double &fp1, const double &fm1, const double &h, const double &eps)
{
    if(fm1==0.0){  cout << " alert!!! " << endl; return fp1/(2.0*h); }  
    
    double dum=fp1/fm1-1.0;
    if(fabs(dum)<=eps) return 0.0; 
    return fm1*dum/(2.0*h);
}

//==================================================================================================
void fcn_for_jacobian_eval(int jX, int neq, double DXj, double z, double Tg, 
                           double NH, Data_Level_I &LI0, double *f)
{
    double Te=Tg;   

    //==============================================================================================
    // reset everything 
    //==============================================================================================
    for(int i=0; i<LI0.neq; i++) f[i]=1.0e-300;        
    
    //==============================================================================================
    LI0.X[jX]+=DXj;
    
    //------------------------------------------------------------------
    // HI & HeI populations
    //------------------------------------------------------------------
    if(jX>=LI0.index_HI && jX<LI0.index_HeI) 
        Hydrogen_Atoms.Set_population_of_level(jX-LI0.index_HI, LI0.X[jX]);
    
    else if(jX>=LI0.index_HeI && jX<LI0.neq-1) 
        HeI_Atoms.Set_Xi(jX-LI0.index_HeI, LI0.X[jX]);       
    
    //------------------------------------------------------------------
    // rho-dependent parts
    //------------------------------------------------------------------
    if(jX==LI0.neq-1)
    {
        double Xp=1.0-Hydrogen_Atoms.X(0), Xe=Xp;
        if(flag_He==1) Xe+=parameters.fHe-HeI_Atoms.Xi(0);
        
        update_Solver_vars_eff(&LI0.X[LI0.index_HI], Hydrogen_Atoms, &LI0.X[LI0.index_HeI], 
                               HeI_Atoms, z, LI0.X[LI0.neq-1], 0);
        
        ODE_HI_effective::evaluate_effective_Rci_Ric_terms_for_Te_derivative(z, Xe, NH*Xp, Ai_df_dx, 
                                                                             Bi_df_dx, get_HI_index, 
                                                                             &f[LI0.index_HI]);
        
        ODE_HI_effective::evaluate_effective_Rij_terms(z, Tg, Hydrogen_Atoms, RijVec_df_dx, 
                                                       get_HI_index, &f[LI0.index_HI]);        
    }
    else 
    {
        //==========================================================================
        // Diffusion correction (added 18.07.2010)
        //==========================================================================
        if(Diffusion_correction_is_on)
        {
            evaluate_HI_Diffusion_correction(z, Hydrogen_Atoms, NH, cosmos.H(z), 
                                             Tg, &f[LI0.index_HI]);
            
            if(DI1_2s_correction_on) 
                evaluate_HI_DI1_DI2_correction(z, Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(1),
                                               Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                               Tg, &f[LI0.index_HI]);
        }      
        
        if(flag_He==1)
        {
            //==================================================================
            // HeI-diffusion correction (added 22.01.2011; JF & JC)
            //==================================================================
            if(Diffusion_correction_HeI_is_on)
            {
                double dum;
                evaluate_HeI_Diffusion_correction(z, HeI_Atoms, NH, cosmos.H(z), Tg, 
                                                  &f[LI0.index_HeI], dum, f[LI0.index_HI]);
                
                if(DI1_2s_correction_on_HeI) 
                    evaluate_HeI_DI1_DI2_correction(z, HeI_Atoms.Xi(0), HeI_Atoms.Xi(1), 
                                                    HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                                    Tg, &f[LI0.index_HeI]);
            }                

            //==========================================================================
            // HeI-feedback
            //==========================================================================
            if(HeISTfeedback==1 && !Diffusion_correction_HeI_is_on)
            {
                if(_HeI_feedback_w_HI_abs==0) 
                    evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                          Dn_HeI_feedback_no_HI_abs, &f[LI0.index_HeI]);
                
                if(_HeI_feedback_w_HI_abs==1) 
                    evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                          Dn_HeI_feedback, &f[LI0.index_HeI]);
            }
            
            //==================================================================
            // absorption by neutral hydrogen
            //==================================================================
            if(flag_HI_absorption==1 && z<=zcrit_HI && !Diffusion_correction_HeI_is_on)
            {  
                if(jX==LI0.index_HeI || jX==LI0.index_HeI+2 || jX==LI0.index_HI)
                {
                    double dXe_HeI_HI_abs_dt=0.0;
                    //==================================================================
                    // HeI Ly-a line
                    //==================================================================
                    evaluate_HI_abs_HeI(z, Hydrogen_Atoms.X(1, 0), NH, cosmos.H(z), Tg, Te,
                                        HeI_Atoms, cosmos, HeI_Atoms.HILyc, 
                                        DP_interpol_S, dXe_HeI_HI_abs_dt);
                    
                    // Helium levels
                    f[LI0.index_HeI+0]+= dXe_HeI_HI_abs_dt;
                    f[LI0.index_HeI+2]+=-dXe_HeI_HI_abs_dt;
                    // HI 1s
                    f[LI0.index_HI]+=-dXe_HeI_HI_abs_dt;
                }
                
                //==================================================================
                // from intercombination line
                //==================================================================
                if(jX==LI0.index_HeI || jX==LI0.index_HI || 
                   jX==LI0.index_HeI+(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1))
                {
                    if(z<=zcrit_HI_Int && flag_spin_forbidden==1)
                    {
                        double dXe_HeI_HI_abs_dt=0.0;
                        evaluate_HI_abs_HeI_Intercombination(z, Hydrogen_Atoms.X(1, 0), NH, 
                                                             cosmos.H(z), Tg, Te, 
                                                             HeI_Atoms, cosmos, HeI_Atoms.HILyc, 
                                                             DP_interpol_T, dXe_HeI_HI_abs_dt);     
                        
                        // Helium levels
                        f[LI0.index_HeI+0]+= dXe_HeI_HI_abs_dt;
                        f[LI0.index_HeI+HeI_Atoms.Get_Level_index(2, 1, 1, 1)]+=-dXe_HeI_HI_abs_dt;
                        // HI 1s
                        f[LI0.index_HI]+=-dXe_HeI_HI_abs_dt;
                    }
                }
            }       
        }
    }
        
    //==============================================================================================
    // DM annihilation
    //==============================================================================================
    if(DM_annihilation && (jX==LI0.index_HI || jX==LI0.index_HeI || jX==LI0.neq-1))
    {
        double Xp=1.0-Hydrogen_Atoms.X(0);
        double XHeII = (flag_He==1 ? parameters.fHe-LI0.X[LI0.index_HeI] : 0.0); 
        
        evaluate_DM_annihilation_terms(z, parameters.fHe, Tg, Xp, XHeII, fDM_CosmoRec, 
                                       f[LI0.index_HI], f[LI0.index_HeI], f[LI0.neq-1]);  
    }           

    //==============================================================================================
    // restore old values
    //==============================================================================================
    LI0.X[jX]-=DXj;
    
    //------------------------------------------------------------------
    // RESET HI & HeI populations
    //------------------------------------------------------------------
    if(jX>=LI0.index_HI && jX<LI0.index_HeI) 
        Hydrogen_Atoms.Set_population_of_level(jX-LI0.index_HI, LI0.X[jX]);
    
    else if(flag_He==1 && jX>=LI0.index_HeI && jX<LI0.neq-1) 
        HeI_Atoms.Set_Xi(jX-LI0.index_HeI, LI0.X[jX]);     
    
    return;
}

//==================================================================================================
void compute_numerical_derivative(int colX, int neq, int neq_HI, int neq_HeI, double z, 
                                  double Tg, double NH, Data_Level_I &LI0, double *r)
{ 
    //===========================================================================
    // r[i] is reset here
    //===========================================================================   
    double X0=LI0.X[colX], DXj=X0*1.0e-7, eps=1.0e-12;
    if(X0==0.0) DXj=1.0e-7;
    
    //===========================================================================
    // numerical derivatives with respect to Xj: here part of the populations and 
    // also recombination rates could be changed
    //===========================================================================
    // get f(Xj+DXj)
    fcn_for_jacobian_eval(colX, neq,  DXj, z, Tg, NH, LI0, Level_III_temp_JAC.g);
    // get f(Xj-DXj)
    fcn_for_jacobian_eval(colX, neq, -DXj, z, Tg, NH, LI0, Level_III_temp_JAC.y);
    
    //===========================================================================
    // derivative with respect to rho --> restore recombination rates after calls
    //===========================================================================
    if(colX==neq-1) 
        update_Solver_vars_eff(&LI0.X[LI0.index_HI], Hydrogen_Atoms, 
                               &LI0.X[LI0.index_HeI], HeI_Atoms, 
                               z, LI0.X[LI0.neq-1], 0);
    
    //===========================================================================
    // define numerical derivative
    //===========================================================================
    r[0]+=df_dx_2point_sym(Level_III_temp_JAC.g[0], Level_III_temp_JAC.y[0], DXj, eps);
    
    r[LI0.index_HI]+=df_dx_2point_sym(Level_III_temp_JAC.g[LI0.index_HI], 
                                      Level_III_temp_JAC.y[LI0.index_HI], 
                                      DXj, eps);
    //
    for(int m=0; m<neq_HI; m++)
    { 
        int kl=LI0.index_HI+get_HI_index(m); 
        r[kl]+=df_dx_2point_sym(Level_III_temp_JAC.g[kl], Level_III_temp_JAC.y[kl], DXj, eps); 
    }
    
    if(flag_He==1)
    {
        r[LI0.index_HeI]+=df_dx_2point_sym(Level_III_temp_JAC.g[LI0.index_HeI], 
                                           Level_III_temp_JAC.y[LI0.index_HeI], 
                                           DXj, eps);
        
        for(int m=0; m<neq_HeI; m++)
        { 
            int kl=LI0.index_HeI+get_HeI_index(m); 
            r[kl]+=df_dx_2point_sym(Level_III_temp_JAC.g[kl], Level_III_temp_JAC.y[kl], DXj, eps); 
        }
    }
    
    return; 
}

//==================================================================================================
//
//
//==================================================================================================
static double jac_effective_z=-1e+300;

//==================================================================================================
//
// jac for ODE evaluation
//
//==================================================================================================
void jac_effective(double z, Data_Level_III &LIII, int col)
{
    int iHI=Level_I_temp.index_HI;
    int iHeI=Level_I_temp.index_HeI;
    int neq_HI=get_number_of_resolved_levels();
    int neq_HeI=get_number_of_resolved_levels_HeI();
    int neq=1+neq_HI+neq_HeI;
    
    //==================================================================
    // first decide if one should update z-dependent functions
    //==================================================================
    if(jac_effective_z!=z)
    {
        jac_effective_z=z;
        //==================================================================
        // transform from yi --> Xi
        //==================================================================
        transform_y2X(LIII, z, Level_I_temp);
        
        //==================================================================
        // copy Xi_H --> Atoms && update recombination coefficients
        //==================================================================
        update_Solver_vars_eff(&Level_I_temp.X[iHI], Hydrogen_Atoms, &Level_I_temp.X[iHeI], 
                               HeI_Atoms, z, Level_I_temp.X[Level_I_temp.neq-1]);
    }   
    
    double rho=Level_I_temp.X[Level_I_temp.neq-1];
    double Tg=cosmos.TCMB(z);
    double Xp=1.0-Level_I_temp.X[iHI];
    double XHeII = (flag_He==1 ? parameters.fHe-Level_I_temp.X[iHeI] : 0.0); 
    double Xe=Xp+XHeII;
    double NH=cosmos.NH(z), H_z=cosmos.H(z), Np=NH*Xp;
    
    //==================================================================
    // reset everything 
    //==================================================================
    for(int i=0; i<Level_I_temp.neq; i++) Level_I_temp.g[i]=LIII.g[i]=0.0;        
    
    //==================================================================
    // Te derivatives
    //==================================================================
    if(col==0)
    {
        ODE_effective::df_drho_evaluate_TM(z, Xe, parameters.fHe, Tg, H_z, 
                                           Level_I_temp.g[Level_I_temp.neq-1]);
        
        compute_numerical_derivative(Level_I_temp.neq-1, neq, neq_HI, neq_HeI, 
                                     z, Tg, NH, Level_I_temp, Level_I_temp.g);
        
        for(int k=0; k<LIII.neq; k++) 
            LIII.g[k]=Level_I_temp.g[k]*f0[Level_I_temp.neq-1]/f0[k];
    }
    //==================================================================
    // HI X1s derivatives
    //==================================================================
    else if(col==1)
    {
        ODE_effective::df_dXHI1s_evaluate_TM(z, Xe, parameters.fHe, rho, Tg, H_z, 
                                             Level_I_temp.g[Level_I_temp.neq-1]);
        
        ODE_effective::df_dX1s_evaluate_2s_two_photon_decay(z, Tg, 
                                                            Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                            const_HI_A2s_1s, Level_I_temp.g[iHI], 
                                                            Level_I_temp.g[iHI+1]); 
        
        ODE_HI_effective::df_dXHI1s_evaluate_effective_Rci_Ric_terms(z, Xe*NH, Np, Ai_df_dx, 
                                                                     Bi_df_dx, get_HI_index, 
                                                                     &Level_I_temp.g[iHI]);
        
        if(flag_He==1) 
            ODE_HeI_effective::df_dXHI1s_evaluate_effective_Rci_Ric_terms(z, XHeII*NH, Ai_df_dx_HeI, 
                                                                          Bi_df_dx_HeI, 
                                                                          get_HeI_index, 
                                                                          &Level_I_temp.g[iHeI]);
        
        int nLy_max=(int)min(Lyn_max_parameter, 
                             Hydrogen_Atoms.Get_n_of_Level(get_HI_index(neq_HI-1)));
        
        for(int n=2; n<=nLy_max; n++)
        {
            Transition_Data_A T=Hydrogen_Atoms.Level(n, 1).Get_Data_A(1, 0);
            
            int nly=Hydrogen_Atoms.Get_Level_index(n, 1);
            ODE_effective::df_dX1s_evaluate_Ly_n_channel(z, Tg, Hydrogen_Atoms.X(0), 
                                                         Hydrogen_Atoms.X(n, 1), NH, H_z, 
                                                         T.A21, T.lambda21, T.Dnu, 
                                                         Level_I_temp.g[iHI], 
                                                         Level_I_temp.g[iHI+nly]);
        }
        
        compute_numerical_derivative(iHI, neq, neq_HI, neq_HeI, z, Tg, NH, 
                                     Level_I_temp, Level_I_temp.g);
        
        for(int k=0; k<LIII.neq; k++) LIII.g[k]=Level_I_temp.g[k]*f0[iHI]/f0[k];
    }
    //==================================================================
    // HI Xi derivatives (i>1s)
    //==================================================================
    else if(col>1 && col<2+neq_HI)
    {
        int m=col-2, ik=get_HI_index(m);
        ///
        if(col==2) 
            ODE_effective::df_dX2s_evaluate_2s_two_photon_decay(z, Tg, 
                                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                                const_HI_A2s_1s, Level_I_temp.g[iHI], 
                                                                Level_I_temp.g[iHI+1]);
        
        ODE_HI_effective::df_dXi_evaluate_effective_Rci_Ric_terms(m, z, Ai_df_dx, Bi_df_dx, 
                                                                  get_HI_index, 
                                                                  &Level_I_temp.g[iHI]);
        
        ODE_HI_effective::df_dXi_evaluate_effective_Rij_terms(m, z, Tg, Hydrogen_Atoms, 
                                                              RijVec_df_dx, get_HI_index, 
                                                              &Level_I_temp.g[iHI]);
        
        int n=Hydrogen_Atoms.Get_n_of_Level(ik);
        int l=Hydrogen_Atoms.Get_l_of_Level(ik);
        
        if(n<=Lyn_max_parameter && l==1)
        {
            Transition_Data_A T=Hydrogen_Atoms.Level(ik).Get_Data_A(1, 0);
            
            ODE_effective::df_dXnp_evaluate_Ly_n_channel(z, Tg, Hydrogen_Atoms.X(0), 
                                                         Hydrogen_Atoms.X(ik), NH, H_z, T.A21, 
                                                         T.lambda21, T.Dnu, Level_I_temp.g[iHI], 
                                                         Level_I_temp.g[iHI+ik]);
            
            //==================================================================
            // nD-1s quadrupole lines (added May 2011)
            //==================================================================
            if(Atom_activate_HI_Quadrupole_lines && n>2 && n<=nD_Quadrupole_max)
                ODE_effective::df_dXnp_evaluate_nD_1s_Q_channel(Hydrogen_Atoms.X(n, 2), 
                                                                Hydrogen_Atoms.Level(n, 2).Get_A21(1, 0), 
                                                                Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 1)], 
                                                                Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 2)]);
        }
        
        //==================================================================
        // nD-1s quadrupole lines (added May 2011)
        //==================================================================
        if(Atom_activate_HI_Quadrupole_lines && n>2 && n<=nD_Quadrupole_max && l==2)
            ODE_effective::df_dXnd_evaluate_nD_1s_Q_channel(Hydrogen_Atoms.X(n, 1), 
                                                            Hydrogen_Atoms.Level(n, 2).Get_A21(1, 0), 
                                                            Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 1)], 
                                                            Level_I_temp.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 2)]);
        
        
        compute_numerical_derivative(iHI+ik, neq, neq_HI, neq_HeI, z, Tg, NH, 
                                     Level_I_temp, Level_I_temp.g);
        
        for(int k=0; k<LIII.neq; k++) LIII.g[k]=Level_I_temp.g[k]*f0[iHI+ik]/f0[k];
    }
    //==================================================================
    // HeI X1s derivatives
    //==================================================================
    else if(flag_He==1 && col==2+neq_HI)
    {
        ODE_effective::df_dXHeI1s_evaluate_TM(z, Xe, parameters.fHe, rho, Tg, H_z, 
                                              Level_I_temp.g[Level_I_temp.neq-1]);
        
        ODE_effective::df_dX1s_evaluate_2s_two_photon_decay(z, Tg, 
                                                            HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                                            const_HeI_A2s_1s, 
                                                            Level_I_temp.g[iHeI], 
                                                            Level_I_temp.g[iHeI+1]); 
        
        ODE_HI_effective ::df_dXHeI1s_evaluate_effective_Rci_Ric_terms(z, Np, Ai_df_dx, Bi_df_dx, 
                                                                       get_HI_index, 
                                                                       &Level_I_temp.g[iHI]);
        
        ODE_HeI_effective::df_dXHeI1s_evaluate_effective_Rci_Ric_terms(z, Xe*NH, XHeII*NH, 
                                                                       Ai_df_dx_HeI, 
                                                                       Bi_df_dx_HeI, 
                                                                       get_HeI_index, 
                                                                       &Level_I_temp.g[iHeI]);
        
        for(int k=1; k<neq_HeI; k++)
        {
            int ik=get_HeI_index(k);
            if(flag_spin_forbidden==0 && ik==(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1)) continue;
            
            Transition_Data_HeI_A T=HeI_Atoms.Get_Trans_Data(ik, 1, 0, 0, 0);
            if(T.A21!=0) 
                ODE_effective::df_dX1s_evaluate_Ly_n_channel(z, Tg, HeI_Atoms.Xi(0), 
                                                             HeI_Atoms.Xi(ik), NH, H_z, T.A21, 
                                                             T.lambda21, T.Dnu, 
                                                             Level_I_temp.g[iHeI], 
                                                             Level_I_temp.g[iHeI+ik], 
                                                             HeI_Atoms.Get_gw(ik)/T.gwp);
        }
        
        compute_numerical_derivative(iHeI, neq, neq_HI, neq_HeI, z, Tg, NH, 
                                     Level_I_temp, Level_I_temp.g);
        
        for(int k=0; k<LIII.neq; k++) LIII.g[k]=Level_I_temp.g[k]*f0[iHeI]/f0[k];
    }
    //==================================================================
    // HeI Xi derivatives (i>1s)
    //==================================================================
    else if(flag_He==1 && col>2+neq_HI)
    {
        int m=col-(3+neq_HI), ik=get_HeI_index(m);
        //
        if(3+neq_HI) 
            ODE_effective::df_dX2s_evaluate_2s_two_photon_decay(z, Tg, 
                                                                HeI_Atoms.Sing.Level(2,0).Get_Dnu_1s2(), 
                                                                const_HeI_A2s_1s, 
                                                                Level_I_temp.g[iHeI], 
                                                                Level_I_temp.g[iHeI+1]); 
        
        ODE_HeI_effective::df_dXi_evaluate_effective_Rci_Ric_terms(m, z, Ai_df_dx_HeI, 
                                                                   Bi_df_dx_HeI, get_HeI_index, 
                                                                   &Level_I_temp.g[iHeI]);
        
        ODE_HeI_effective::df_dXi_evaluate_effective_Rij_terms(m, z, Tg, HeI_Atoms, 
                                                               RijVec_df_dx_HeI, 
                                                               get_HeI_index, 
                                                               &Level_I_temp.g[iHeI]);        
        
        int l=HeI_Atoms.Get_l(ik);
        int J=HeI_Atoms.Get_J(ik);
        if((l==1 && J==1) || (l==2 && J==2))
        {
            if(!( flag_spin_forbidden==0 && ik==(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1) ) )
            {
                Transition_Data_HeI_A T=HeI_Atoms.Get_Trans_Data(ik, 1, 0, 0, 0);
                
                ODE_effective::df_dXnp_evaluate_Ly_n_channel(z, Tg, HeI_Atoms.Xi(0), 
                                                             HeI_Atoms.Xi(ik), NH, H_z, 
                                                             T.A21, T.lambda21, T.Dnu, 
                                                             Level_I_temp.g[iHeI], 
                                                             Level_I_temp.g[iHeI+ik], 
                                                             HeI_Atoms.Get_gw(ik)/T.gwp);
            }
        }       
        
        compute_numerical_derivative(iHeI+ik, neq, neq_HI, neq_HeI, z, Tg, NH, 
                                     Level_I_temp, Level_I_temp.g);
        
        for(int k=0; k<LIII.neq; k++) LIII.g[k]=Level_I_temp.g[k]*f0[iHeI+ik]/f0[k];
    }
    else if(flag_He==0 && col>=2+neq_HI){}
    else{ cerr << " This should not happen. Exiting. " << endl; exit(0); }
    
    return;
}

//==================================================================================================
void fcn_effective(int *neq, double *z, double *y, double *f, int col)
{    
    int index_HI=Level_III_temp.index_HI;
    int index_HeI=Level_III_temp.index_HeI;
    int nres_HI=get_number_of_resolved_levels();
    int nres_HeI=get_number_of_resolved_levels_HeI();
        
    //==================================================================
    copy_ysol_to_LIII(y, Level_III_temp.y, Level_III_temp.neq, index_HI, index_HeI);
    
    //==================================================================
    // col<0 : evaluate the rhs of system  dy/dz == g(t,y); 
    // col>=0: evaluate derivative with respect to j==col
    //==================================================================
    if(col<0) fcn_effective(*z, Level_III_temp); 
    else jac_effective(*z, Level_III_temp, col); 
    
    //==================================================================
    // remember that Level_III_temp.g is g(t, y) --> need dz/dt factor!!!
    //==================================================================
    double dz_dt=-cosmos.H( (*z) )*( 1.0+ (*z) );
    
    //==================================================================
    f[0]=Level_III_temp.g[Level_III_temp.neq-1]/dz_dt;
    f[1]=Level_III_temp.g[index_HI]/dz_dt;
    
    for(int k=0; k<nres_HI; k++) f[k+2]=Level_III_temp.g[index_HI+get_HI_index(k)]/dz_dt;
    //
    if(flag_He==1)
    {
        f[nres_HI+2]=Level_III_temp.g[index_HeI]/dz_dt;
        
        for(int k=0; k<nres_HeI; k++) 
            f[nres_HI+k+3]=Level_III_temp.g[index_HeI+get_HeI_index(k)]/dz_dt;
    }
    
    return;
}   

//==================================================================================================
//==================================================================================================
