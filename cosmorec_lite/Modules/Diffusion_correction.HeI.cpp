//====================================================================================================================
//
// for module to account for the diffusion corrections to escape probabilities
//
//====================================================================================================================
DPesc_splines_Data DPesc_splines_Data_HeI[50];
int DPesc_splines_HeI_are_setup=0;

vector<int> DPesc_HeI_level_indices;


//====================================================================================================================
//
// for module to account for the HeI 2s-1s corrections
//
//====================================================================================================================
bool DI1_2s_correction_on_HeI=0;


//====================================================================================================================
//
// This part is just to set up interpolations for the computed correction functions from the 
// diffusion code. Splines are used currently, but linear interpolation might work better 
// when sampling in redshift is sparser.
//
//====================================================================================================================
void setup_DP_spline_data_HeI(int n, vector<double> &DPesc_vec_z, vector<double> &DPesc_vec)
{
    int nz=DPesc_vec_z.size();
    //
    double *za=new double[nz];
    double *ya=new double[nz];
    
    for(int i=0; i<nz; i++)
    { 
        za[i]=DPesc_vec_z[nz-1-i]; 
        ya[i]=DPesc_vec[nz-1-i]; 
    }
    
    if(DPesc_splines_HeI_are_setup==0) 
        DPesc_splines_Data_HeI[n].memindex=calc_spline_coeffies_JC(nz, za, ya,
                                                                   "setup_DP_spline_data_HeI :: "
                                                                   +int_to_string(n));   
    //
    else update_spline_coeffies_JC(DPesc_splines_Data_HeI[n].memindex, nz, za, ya);                

    DPesc_splines_Data_HeI[n].zmin=max(1300.0, DPesc_vec_z[nz-1]);
    DPesc_splines_Data_HeI[n].zmax=min(3200.0, DPesc_vec_z[0]);
    DPesc_splines_Data_HeI[n].DPmin=DPesc_splines_Data_HeI[n].DPmax=0.0;
    
    delete [] za;
    delete [] ya;
    
    return;
}

//====================================================================================================================
void setup_DP_spline_data_HeI(vector<double> &DF_vec_z, 
                              vector<vector<double> > &DF_vec, 
                              vector<double> &DI1_2s_vec,
                              vector<int> &HeI_level_indices)
{   
    if(DI1_2s_correction_on_HeI) setup_DP_spline_data_HeI(49, DF_vec_z, DI1_2s_vec);

    DPesc_HeI_level_indices=HeI_level_indices;
    
    for(int m=0; m<(int)DF_vec.size(); m++) setup_DP_spline_data_HeI(m, DF_vec_z, DF_vec[m]);

    DPesc_splines_HeI_are_setup=1;

    return;
}

//====================================================================================================================
//
// Setup interpolation functions for rate corrections
//
//====================================================================================================================
void setup_DF_interpol_data_HeI(vector<double> &DF_vec_z, 
                                vector<vector<double> > &DF_vec,
                                vector<double> &DI1_2s_vec,
                                vector<int> &HeI_level_indices)
{
    setup_DP_spline_data_HeI(DF_vec_z, DF_vec, DI1_2s_vec, HeI_level_indices);
    
    return;
}

//====================================================================================================================
//
// mapping of functions
//
//====================================================================================================================
double interpolate_DF_HeI(double z, int k)
{
    if(z<=DPesc_splines_Data_HeI[k].zmin) return DPesc_splines_Data_HeI[k].DPmin;
    if(z>=DPesc_splines_Data_HeI[k].zmax) return DPesc_splines_Data_HeI[k].DPmax;
    
    return calc_spline_JC(z, DPesc_splines_Data_HeI[k].memindex);
}

//====================================================================================================================
double DP_spline_HeI(int i, double z){ return interpolate_DF_HeI(z, i); }
double DI1_2s_spline_HeI(double z){ return interpolate_DF_HeI(z, 49); }

//====================================================================================================================
//
// correction to HeI Singlet 2s-1s two-photon channel
//
//====================================================================================================================
void evaluate_HeI_DI1_DI2_correction(double z, double X1s, double X2s, double nu21, double Tg, double *g)
{
    double x21=const_h_kb*nu21/Tg;  
    double DI1_2s=(DI1_2s_correction_on_HeI==1 ? DI1_2s_spline_HeI(z) : 0.0)*(X2s-X1s*exp(-x21));
    double DRij=const_HeI_A2s_1s*DI1_2s;
    
    g[0]+=-DRij; // remove electron from 1s level
    g[1]+= DRij; // add it to 2s-state
    
//  cout << z << " HeI 2s-1s correction " << DRij << " DI1= " << DI1_2s << " " << exp(-x21) << endl;
    
    return;
}

//====================================================================================================================
//
// add HeI Diffusion correction
//
//====================================================================================================================
void evaluate_HeI_Diffusion_correction(double z, Gas_of_HeI_Atoms &HeIA, 
                                       double NH, double H_z, double Tg, 
                                       double *g_HeI, double &g_Xe, double &g_HI1s)
{
//    return;
    
    double DRij;
    double X1s=HeIA.Xi(0);
    int i_HeI;
    Transition_Data_HeI_A T;
    
    for(int line=0; line<(int)DPesc_HeI_level_indices.size(); line++)
    {
        i_HeI=DPesc_HeI_level_indices[line];
        T=HeIA.Get_Trans_Data(i_HeI, 1, 0, 0, 0);
        double w=HeIA.Get_gw(i_HeI)/T.gwp;
        double exp_x=exp(-const_h_kb*T.Dnu/Tg);
        double DPval=DP_spline_HeI(line, z);

        DRij=T.A21*1.0/(1.0-exp_x)*DPval*w*X1s*(HeIA.Xi(i_HeI)/X1s/w-exp_x);

        //================================================================================
        // Helium transitions
        //================================================================================
        g_HeI[0]      += DRij;   // add electron to HeI 1s state
        g_HeI[i_HeI]  +=-DRij;   // remove it from upper HeI state
        //================================================================================

        //================================================================================
        // corresponding transitions between continuum and HI ground state
        //================================================================================
        g_Xe          += DRij;   // add electron to continuum (HI ionization)
        g_HI1s        +=-DRij;   // remove it from HI 1s state       
        //================================================================================
    }
    
    return;
}

//====================================================================================================================
//====================================================================================================================
