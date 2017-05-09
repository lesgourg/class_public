//==================================================================================================
//
// Module to account for the diffusion corrections to escape probabilities
//
//==================================================================================================
int DPesc_splines_are_setup=0;

struct DPesc_splines_Data
{
    int memindex;
    double zmin, zmax;
    double DPmin, DPmax;
    vector<double> za, ya; 
    int accel;
};

DPesc_splines_Data DPesc_splines_Data_DF[100];

//==================================================================================================
int index_2g_spline;
int index_R_spline;
int index_Ly_n_spline;
int index_nD1s_spline;
int index_all_spline;

#define DF_SPLINES

//#define DIFF_CORR_STOREII   // Changes between two ways of storing the diffusion correction.
                            // Not setting it is consistent with CosmoRec v1.5
                            // make sure to also change this in ./PDE_Problem/Solve_PDEs_integrals.cpp

//==================================================================================================
//
// for module to account for the 2s-1s corrections
//
//==================================================================================================
int induced_flag=0;          // is set according to parameter file (see ./Modules/aux_functions.cpp)
bool DI1_2s_correction_on=0;


//==================================================================================================
// 
// This part is just to set up interpolations for the computed correction functions from the 
// diffusion code. Splines are used currently, but linear interpolation might work better 
// when sampling in redshift is sparser.
//
//==================================================================================================
void setup_DF_spline_data(int n, vector<double> &DPesc_vec_z, vector<double> &DPesc_vec)
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
    
    if(DPesc_splines_are_setup==0) 
        DPesc_splines_Data_DF[n].memindex=calc_spline_coeffies_JC(nz, za, ya,
                                                                  "setup_DF_spline_data :: "
                                                                  +int_to_string(n));   
    //
    else update_spline_coeffies_JC(DPesc_splines_Data_DF[n].memindex, nz, za, ya);                

    DPesc_splines_Data_DF[n].zmin=max(Diff_corr_zmin, DPesc_vec_z[nz-1]);
    DPesc_splines_Data_DF[n].zmax=min(Diff_corr_zmax, DPesc_vec_z[0]);
    DPesc_splines_Data_DF[n].DPmin=DPesc_splines_Data_DF[n].DPmax=0.0;
    
    delete [] za;
    delete [] ya;
    
    return;
}

//==================================================================================================
void setup_DF_spline_data(vector<double> &DF_vec_z, 
                          vector<vector<double> > &DF_Ly_n_vec, 
                          vector<vector<double> > &DF_2g_vec, 
                          vector<vector<double> > &DF_R_vec, 
                          vector<vector<double> > &DF_nD1s_vec, 
                          vector<double> &DI1_2s_vec)
{   
    if(DI1_2s_correction_on) setup_DF_spline_data(0, DF_vec_z, DI1_2s_vec);

    index_2g_spline=2;
    //
    for(int m=0; m<(int)DF_2g_vec.size(); m++) 
        setup_DF_spline_data(m+index_2g_spline, DF_vec_z, DF_2g_vec[m]);
    index_R_spline=index_2g_spline+DF_2g_vec.size();
    //
    for(int m=0; m<(int)DF_R_vec.size(); m++) 
        setup_DF_spline_data(m+index_R_spline, DF_vec_z, DF_R_vec[m]);
    index_Ly_n_spline=index_R_spline+DF_R_vec.size();
    //
    for(int m=0; m<(int)DF_Ly_n_vec.size(); m++) 
        setup_DF_spline_data(m+index_Ly_n_spline, DF_vec_z, DF_Ly_n_vec[m]);
    index_nD1s_spline=index_Ly_n_spline+DF_Ly_n_vec.size();
    //
    for(int m=0; m<(int)DF_nD1s_vec.size(); m++) 
        setup_DF_spline_data(m+index_nD1s_spline, DF_vec_z, DF_nD1s_vec[m]);
    index_all_spline=index_nD1s_spline+DF_nD1s_vec.size();
    //
    DPesc_splines_are_setup=1;

    return;
}

//==================================================================================================
//
// n-point interpolation option
//
//==================================================================================================
void setup_DF_n_point_data(int n, vector<double> &DPesc_vec_z, vector<double> &DPesc_vec)
{
    int nz=DPesc_vec_z.size();
    //
    DPesc_splines_Data_DF[n].za=DPesc_vec_z;
    DPesc_splines_Data_DF[n].ya=DPesc_vec;
    DPesc_splines_Data_DF[n].memindex=n;
    DPesc_splines_Data_DF[n].zmin=max(Diff_corr_zmin, DPesc_vec_z[nz-1]);
    DPesc_splines_Data_DF[n].zmax=min(Diff_corr_zmax, DPesc_vec_z[0]);
    DPesc_splines_Data_DF[n].DPmin=DPesc_splines_Data_DF[n].DPmax=0.0;
    DPesc_splines_Data_DF[n].accel=0;
    
    return;
}

//==================================================================================================
void setup_DF_n_point_data(vector<double> &DF_vec_z, 
                           vector<vector<double> > &DF_Ly_n_vec, 
                           vector<vector<double> > &DF_2g_vec, 
                           vector<vector<double> > &DF_R_vec, 
                           vector<vector<double> > &DF_nD1s_vec, 
                           vector<double> &DI1_2s_vec)
{   
    if(DI1_2s_correction_on) setup_DF_n_point_data(0, DF_vec_z, DI1_2s_vec);
    
    index_2g_spline=2;
    //
    for(int m=0; m<(int)DF_2g_vec.size(); m++) 
        setup_DF_n_point_data(m+index_2g_spline, DF_vec_z, DF_2g_vec[m]);
    index_R_spline=index_2g_spline+DF_2g_vec.size();
    //
    for(int m=0; m<(int)DF_R_vec.size(); m++) 
        setup_DF_n_point_data(m+index_R_spline, DF_vec_z, DF_R_vec[m]);
    index_Ly_n_spline=index_R_spline+DF_R_vec.size();
    //
    for(int m=0; m<(int)DF_Ly_n_vec.size(); m++) 
        setup_DF_n_point_data(m+index_Ly_n_spline, DF_vec_z, DF_Ly_n_vec[m]);
    index_nD1s_spline=index_Ly_n_spline+DF_Ly_n_vec.size();
    //
    for(int m=0; m<(int)DF_nD1s_vec.size(); m++) 
        setup_DF_n_point_data(m+index_nD1s_spline, DF_vec_z, DF_nD1s_vec[m]);
    index_all_spline=index_nD1s_spline+DF_nD1s_vec.size();
    //
    DPesc_splines_are_setup=1;
    
    return;
}

//==================================================================================================
//
// Setup interpolation functions for rate corrections
//
//==================================================================================================
void setup_DF_interpol_data(vector<double> &DF_vec_z, 
                            vector<vector<double> > &DF_Ly_n_vec, 
                            vector<vector<double> > &DF_2g_vec, 
                            vector<vector<double> > &DF_R_vec, 
                            vector<vector<double> > &DF_nD1s_vec, 
                            vector<double> &DI1_2s_vec)
{   
#ifdef DF_SPLINES
    setup_DF_spline_data(DF_vec_z, DF_Ly_n_vec, DF_2g_vec, DF_R_vec, DF_nD1s_vec, DI1_2s_vec);
#else
    setup_DF_n_point_data(DF_vec_z, DF_Ly_n_vec, DF_2g_vec, DF_R_vec, DF_nD1s_vec, DI1_2s_vec);
#endif
    return;
}

//==================================================================================================
//
// mapping of functions
//
//==================================================================================================
double interpolate_DF(double z, int k)
{
    if(z<=DPesc_splines_Data_DF[k].zmin) return DPesc_splines_Data_DF[k].DPmin;
    if(z>=DPesc_splines_Data_DF[k].zmax) return DPesc_splines_Data_DF[k].DPmax;

#ifdef DF_SPLINES
    return calc_spline_JC(z, DPesc_splines_Data_DF[k].memindex);
#else
    double y, dy;
    polint_JC(&DPesc_splines_Data_DF[k].za[0], &DPesc_splines_Data_DF[k].ya[0], 
              DPesc_splines_Data_DF[k].za.size(), z, DPesc_splines_Data_DF[k].accel, 
              2, &y, &dy);
    return y;
#endif
}

//==================================================================================================
double DF_2g_spline(int i, double z){ return interpolate_DF(z, index_2g_spline+i); }
double DF_R_spline(int i, double z){ return interpolate_DF(z, index_R_spline+i); }
double DF_Ly_n_spline(int i, double z){ return interpolate_DF(z, index_Ly_n_spline+i); }
double DF_nD1s_spline(int i, double z){ return interpolate_DF(z, index_nD1s_spline+i); }
double DI1_2s_spline(double z){ return interpolate_DF(z, 0); }

//==================================================================================================
//
// correction to 2s-1s two-photon channel
//
//==================================================================================================
void evaluate_HI_DI1_DI2_correction(double z, double X1s, double X2s, double nu21, 
                                    double Tg, double *g)
{
#ifdef DIFF_CORR_STOREII
    double x21=const_h_kb*nu21/Tg;
    double DI1_2s=DI1_2s_spline(z)*(X2s-X1s*exp(-x21));
    double DRij=const_HI_A2s_1s*DI1_2s;
#else
    double DI1_2s=DI1_2s_spline(z);
    double DRij=DI1_2s*X1s;
#endif
    
    g[0]+=-DRij; // remove electron from 1s level
    g[1]+= DRij; // add it to 2s-state
    
//    cout << z << " 2s-1s correction " << DRij << " DI1= " << DI1_2s << " " << exp(-x21) << endl;
    
    return;
}

//==================================================================================================
//
// Ly-n 1+1 photon correction
//
//==================================================================================================
double compute_DR_Ly_n_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni)
{
    double DRij=DF_Ly_n_spline(i, z)*( HIA.X(ni, 1)/X1s/3.0-exp(-const_h_kb*HIA.Level(ni, 1).Get_Dnu_1s()/Tg) );

//    cout << z << " (ni, li) " << ni << " " << 1 << " Ly-n " << DRij << endl;
    return DRij;
}

void evaluate_HI_Diffusion_correction_Ly_n(double z, Gas_of_Atoms &HIA, 
                                           double NH, double H_z, double Tg, double *g)
{
    if(index_Ly_n_spline==index_nD1s_spline) return;
    
    double DRij;
    double X1s=HIA.X(1, 0);
    int index=1;
    int nmax=8;
    
    for(int ni=3; ni<=nmax && index+index_Ly_n_spline<index_nD1s_spline; ni++)
    {
        DRij=compute_DR_Ly_n_func(index, z, X1s, Tg, HIA, ni);
        g[0]+=-DRij;                          // remove electron from lower level
        g[HIA.Get_Level_index(ni, 1)]+= DRij; // add it to -state
        index++;
    }
    
    return;
}

//==================================================================================================
//
// two-photon correction
//
//==================================================================================================
double compute_DR_2g_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni, int li)
{
#ifdef DIFF_CORR_STOREII
    double DRij=DF_2g_spline(i, z)*( HIA.X(ni, li)/X1s/(2.0*li+1.0)-exp(-const_h_kb*HIA.Level(ni, li).Get_Dnu_1s()/Tg) );
#else
    double DRij=X1s*DF_2g_spline(i, z);
#endif
    
//    cout << z << " (ni, li) " << ni << " " << li << " 2g " << DRij << endl;
    return DRij;
}

void evaluate_HI_Diffusion_correction_2_gamma(double z, Gas_of_Atoms &HIA, 
                                              double NH, double H_z, double Tg, double *g)
{
    if(index_2g_spline==index_R_spline) return;

    double DRij;
    double X1s=HIA.X(1, 0);
    int index=0;
    int nmax=8;
    int lmax=2;

    for(int ni=3; ni<=nmax && index+index_2g_spline<index_R_spline; ni++)
        for(int li=0; li<=lmax; li+=2)
        {
            DRij=compute_DR_2g_func(index, z, X1s, Tg, HIA, ni, li);
            g[0]+=-DRij;                           // remove electron from lower level
            g[HIA.Get_Level_index(ni, li)]+= DRij; // add it to -state
            index++;
        }
    
    return;
}

//==================================================================================================
//
// Raman correction
//
//==================================================================================================
double compute_DR_R_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni, int li)
{
#ifdef DIFF_CORR_STOREII
    double DRij=DF_R_spline(i, z)*( HIA.X(ni, li)/X1s/(2.0*li+1.0)-exp(-const_h_kb*HIA.Level(ni, li).Get_Dnu_1s()/Tg) );
#else
    double DRij=X1s*DF_R_spline(i, z);
#endif

//    cout << z << " (ni, li) " << ni << " " << li << " Raman " << DRij << endl;
    return DRij;
}

void evaluate_HI_Diffusion_correction_R(double z, Gas_of_Atoms &HIA, 
                                        double NH, double H_z, double Tg, double *g)
{
    if(index_R_spline==index_Ly_n_spline) return;
    
    double DRij;
    double X1s=HIA.X(1, 0);
    int index=0;
    int nmax=7;
    int lmax=2;
    
    for(int ni=2; ni<=nmax && index+index_R_spline<index_Ly_n_spline; ni++)
        for(int li=0; li<=lmax && li<ni; li+=2)
        {
            DRij=compute_DR_R_func(index, z, X1s, Tg, HIA, ni, li);
            g[0]+=-DRij;                           // remove electron from lower level
            g[HIA.Get_Level_index(ni, li)]+= DRij; // add it to -state
            index++;
        }
    
    return;
}

//==================================================================================================
// 
// nD-1s quadrupole lines 1+1 photon correction (added June 2011)
//
//==================================================================================================
double compute_DR_nD1s_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni)
{
    double DRij=DF_nD1s_spline(i, z)*( HIA.X(ni, 2)/X1s/5.0-exp(-const_h_kb*HIA.Level(ni, 1).Get_Dnu_1s()/Tg) );

//    cout << z << " (ni, li) " << ni << " " << 2 << " nD-1s " << DRij << endl;
    return DRij;
}

void evaluate_HI_Diffusion_correction_nD1s(double z, Gas_of_Atoms &HIA, 
                                           double NH, double H_z, double Tg, double *g)
{
    if(index_nD1s_spline==index_all_spline) return;
    
    double DRij;
    double X1s=HIA.X(1, 0);
    int index=0;
    int nmax=8;
    
    for(int ni=3; ni<=nmax && index+index_nD1s_spline<index_all_spline; ni++)
    {
        DRij=compute_DR_nD1s_func(index, z, X1s, Tg, HIA, ni);
        g[0]+=-DRij;                          // remove electron from lower level
        g[HIA.Get_Level_index(ni, 2)]+= DRij; // add it to -state
        index++;
    }
    
    return;
}

//==================================================================================================
//
// change this part to switch play with the different correction terms; if the Ly-n correction is 
// included the Raman and two-photon correction should not be activated.
// 
//==================================================================================================
void evaluate_HI_Diffusion_correction(double z, Gas_of_Atoms &HIA, 
                                      double NH, double H_z, double Tg, double *g)
{
    evaluate_HI_Diffusion_correction_2_gamma(z, HIA, NH, H_z, Tg, g);
    evaluate_HI_Diffusion_correction_R      (z, HIA, NH, H_z, Tg, g);
//    evaluate_HI_Diffusion_correction_Ly_n   (z, HIA, NH, H_z, Tg, g);
//    evaluate_HI_Diffusion_correction_nD1s   (z, HIA, NH, H_z, Tg, g);
    
    return;
}

//==================================================================================================
//==================================================================================================
