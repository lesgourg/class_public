//===========================================================================================================
// some flags for the solver
//===========================================================================================================
int rescale=1;

//===========================================================================================================
//
// Level I <--> III
//
//===========================================================================================================

//===========================================================================================================
// variable: ( si=Xi/fi )
//===========================================================================================================
double X2y(int kX, const Data_Level_I &LI, double z)
{ return LI.X[kX]/f0[kX]; }

double y2X(int kX, const Data_Level_III &LIII, double z)
{  return LIII.y[kX]*f0[kX]; }

double g_tX2g_ty(int kX, const Data_Level_I &LI, double z){ return LI.g[kX]/f0[kX]; }

//===========================================================================================================
void transform_y2X(const Data_Level_III &LIII, double z, Data_Level_I &LI)
{ 
  for(int k=0; k<LI.neq; k++) LI.X[k]=y2X(k, LIII, z);
  return;
}

//===========================================================================================================
void transform_g_tX2g_ty(const Data_Level_I &LI, double z, Data_Level_III &LIII)
{ 
  for(int k=0; k<LIII.neq; k++) LIII.g[k]=g_tX2g_ty(k, LI, z);
  return;
}

//===========================================================================================================
//
// initialize conversion factors
//
//===========================================================================================================
void init_const_factors_1(double *ff0, int neq)
{ 
    for(int k=0; k<neq; k++) ff0[k]=1.0;
    
    return;
}

//===========================================================================================================
// for reference levels
//===========================================================================================================
const double const_kappa_H_var_trans=const_EH_inf_ergs/(1.0+const_me_mp)/const_kB; 

double DE_ij_kbTg(int ni, int nj, double z)                   // ni > nj >=1
{ return const_kappa_H_var_trans*(1.0/nj/nj-1.0/ni/ni)/cosmos.TCMB(z); }

double DE_kbTg_HeI(double Dnu, double z)  
{ return const_h*Dnu/const_kB/cosmos.TCMB(z); }

//===========================================================================================================
void init_const_factors_2(double *ff0, const Data_Level_III &LIII, 
                          Gas_of_Atoms &H_Atoms, Gas_of_HeI_Atoms &HeIA, 
                          double zH=zref, double zHe=zref_HeI)
{
    //=================================================================================
    // Xi=f0i*si
    //=================================================================================
    ff0[0]=1.0;                             // Xe
    ff0[1]=1.0;                             // X1s
    ff0[2]=1.0e-14;                         // X2s
    ff0[3]=3.0e-14;                         // X2p
    ff0[Level_I.neq-1]=1.0;                 // rho  
    ff0[1]=1.0e-4;                          // X1s

    //=================================================================================
    // all the other levels are rescaled using the Saha-exponential relative to the 2s
    //=================================================================================
    for(int k=4; k<Level_I.index_HeI; k++) 
    {
        ff0[k]=gi_f(H_Atoms.Get_l_of_Level(k-1))/gi_f(0)*exp(-DE_ij_kbTg(H_Atoms.Get_n_of_Level(k-1), 2, zH) )*ff0[2];
    }
    
    //=================================================================================
    // helium part
    //=================================================================================
    int iTrip=HeIA.Get_indexT(), iTripnoj=HeIA.Get_indexT_no_j();
    double nuion_2s_S=HeIA.Sing.Level(2, 0).Get_nu_ion();
    double nuion_2s_T=HeIA.Trip.Level(2, 0, 1).Get_nu_ion();
    double nuion;
    
    //=================================================================================
    // Singlet 1s, 2s, 2p
    //=================================================================================
    ff0[LIII.index_HeI+0]=1.0;
    ff0[LIII.index_HeI+1]=1.0e-16;
    ff0[LIII.index_HeI+2]=3.0e-16;
    for(int k=3; k<(int)HeIA.Sing.Get_total_number_of_Levels(); k++)
    {
        nuion=HeIA.Sing.Level(k).Get_nu_ion();
        ff0[k+LIII.index_HeI]=HeIA.Get_gw(k)/1.0
        *exp( -DE_kbTg_HeI(fabs(nuion_2s_S-nuion), zHe) )*ff0[LIII.index_HeI+1];
    }
    
    //=================================================================================
    // Triplet 2s, 2p (j-resolved)
    //=================================================================================
    ff0[LIII.index_HeI+iTrip+0]=3.0e-15;
    ff0[LIII.index_HeI+iTrip+1]=1.0*ff0[LIII.index_HeI+iTrip+0];
    ff0[LIII.index_HeI+iTrip+2]=3.0*ff0[LIII.index_HeI+iTrip+0];
    ff0[LIII.index_HeI+iTrip+3]=5.0*ff0[LIII.index_HeI+iTrip+0];
    
    for(int k=4; k<(int)HeIA.Trip.Get_total_number_of_Levels(); k++)
    {
        nuion=HeIA.Trip.Level(k).Get_nu_ion();
        ff0[LIII.index_HeI+iTrip+k]=HeIA.Get_gw(k+iTrip)/3.0
        *exp( -DE_kbTg_HeI(fabs(nuion_2s_T-nuion), zHe) )*ff0[LIII.index_HeI+iTrip];      
    }
    
    //=================================================================================
    // Triplet (non-j-resolved)
    //=================================================================================
    for(int k=0; k<(int)HeIA.Trip_no_j.Get_total_number_of_Levels(); k++)
    {
        nuion=HeIA.Trip_no_j.Level(k).Get_nu_ion();
        ff0[LIII.index_HeI+iTripnoj+k]=HeIA.Get_gw(k+iTripnoj)/3.0
        *exp( -DE_kbTg_HeI(fabs(nuion_2s_T-nuion), zHe) )*ff0[LIII.index_HeI+iTrip];      
    }
    
    return;
}

//===========================================================================================================
void set_variables(double zH=zref, double zHe=zref_HeI)
{
    if(rescale==0) init_const_factors_1(f0, Level_III.neq);
    else if(rescale==1) init_const_factors_2(f0, Level_III, Hydrogen_Atoms, HeI_Atoms, zH, zHe);

    return;
}

//===========================================================================================================
//
// some auxillary routines for all variables 
// below the Data structures are not global!!!
//
//===========================================================================================================
void copy_y_to_X(const Data_Level_III &LIII, double z, Data_Level_I &LI)
{
    for(int k=0; k<LI.neq; k++) LI.X[k]= y2X(k, LIII, z);   
    return;
}

//===========================================================================================================
void initialize_y(const Data_Level_I &LI, Gas_of_Atoms &H_Atoms, Gas_of_HeI_Atoms &HeIA, 
                  double z, Data_Level_III &LIII, int wait=0)
{
    int mess=0;
    int iTrip=HeIA.Get_indexT(), iTripnoj=HeIA.Get_indexT_no_j();
    for(int k=0; k<LIII.neq; k++) 
    {
        if(mess==1)
        {
            if(k==LIII.index_HI) cout << " Hydrogen levels " << endl;
            if(k==LIII.index_HeI && flag_He==1) cout << " Helium levels (Singlet) " << endl;
            if(k==LIII.index_HeI+iTrip && flag_He==1) cout << " Helium levels (Triplet j-resolved) " << endl;
            if(k==LIII.index_HeI+iTripnoj && njresolved<nShellsHeI && flag_He==1) cout << " Helium levels (Triplet non-j-resolved)" << endl;
        }
        
        if(flag_He==0 && LIII.index_HeI<=k && k<LIII.index_HeI+LIII.nHeIeq) continue;
        
        LIII.y[k]=X2y(k, LI, z);
        if(mess==1) cout << " initialize_y: " << k << " " << LIII.y[k] << " " << f0[k] << endl;
    }
    
    if(wait!=0) wait_f_r();
    
    return;
}

//===========================================================================================================
void show_results_y(double z, const Data_Level_III &LIII)
{
    cout << " z: " << z << " ye: " << LIII.y[0] << " rho: "<< LIII.y[LIII.neq-1] << " y_i: ";
    for(int k=1; k<LIII.neq-1; k++) cout << LIII.y[k] << " ";
    cout << endl;
    
    return;
}

//===========================================================================================================
void show_results_X(double z, const Data_Level_I &LI)
{
    cout << " z: " << z << " Xe: " << LI.X[0] << " rho: " << LI.X[LI.neq-1] << " X_i: ";
    for(int k=1; k<LI.neq-1; k++) cout << LI.X[k] << " ";
    
    cout << endl;
    
    return;
}

