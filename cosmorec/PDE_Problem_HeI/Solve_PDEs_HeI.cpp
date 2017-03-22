//====================================================================================================================
// Authors Jens Chluba & Jeffrey Fung Feb-May 2011
//====================================================================================================================
// 30.07.2014: Xi_Data is passed on read only

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <sstream>

//====================================================================================================================
// libs and routines
//====================================================================================================================
#include "physical_consts.h"
#include "routines.h"
#include "Cosmos.h"

#include "Sobolev.h"
#include "Patterson.h"

#include "get_effective_rates.HeI.h"
#include "Load.Populations.HI.h"
#include "Load.Populations.HeI.h"
#include "HeI_pd_Rp_splines_effective.h"

#include "PDE_solver.h"
#include "define_PDE_HeI.h"
#include "Solve_PDEs_HeI.h"

//====================================================================================================================
using namespace std;

string HeI_outputdir=COSMORECDIR+"./temp/";

//====================================================================================================================
//
// flags for message output and aVoigt setting
//
//====================================================================================================================
int show_messages_HeI=-1;


//====================================================================================================================
//
// settings for the grid
//
//====================================================================================================================
int core_HeI=25;            // points per Doppler width in the core region; 
                            // For +/- 10 xD core=25 implies 500 points (which is recommended).
int dec_HeI=50;             // points per decade; from 10^-4 --> 1 for ndec==50 gives ~ 200 points

double theta_HeI=0.55;


//====================================================================================================================
//
// switch on/off electron scattering kernel corrections
//
//====================================================================================================================
int solve_with_Kernel=0; // ==0: no kernel correction; 
                         // ==1: kernel correction using iteration; This runmode takes rather long, but the additional
                         //      correction is small. Still with openmp support it may be reasonable to use in the 
                         //      final iteration of the PDE solver.

//====================================================================================================================
void switch_off_e_scat_kernel_HeI(){ solve_with_Kernel=0; return; }
void switch_on_e_scat_kernel_HeI (){ solve_with_Kernel=1; return; }


//====================================================================================================================
//
// setting different paramters for PDE problem regarding the HeI Singlet 2s-1s channel
//
//====================================================================================================================
bool compute_full_2s_1s_corr_HeI=1;
bool HeI_2s_1s_corr_is_on=1;

//====================================================================================================================
// switch on/off HeI singlet 2s-1s channel; this overwrites all other settings for 2s-1s channel
//====================================================================================================================
int npts_2s1s_HeI=500;        // number of points for HeI singlet 2s-1s two-photon part
double xmin_glob_HeI=0.005;   // if 2s-1s correction is on, this should be <~0.3

//====================================================================================================================
void switch_off_1s_2s_HeI_correction()
{ 
    HeI_switch_off_2s1s(); 
    HeI_2s_1s_corr_is_on=0; 
    
    npts_2s1s_HeI=250;
    xmin_glob_HeI=0.8;
    
    return; 
}

void switch_on_1s_2s_HeI_correction()
{ 
    HeI_switch_on_2s1s(); 
    HeI_2s_1s_corr_is_on=1; 
    
    npts_2s1s_HeI=500;
    xmin_glob_HeI=0.005;
    
    return; 
}

//====================================================================================================================
// HI 1s-2s absorption can be switched off separately with this.
//====================================================================================================================
void switch_off_1s_2s_HeI_absorption_correction()
{
    compute_full_2s_1s_corr_HeI=0;
    return;
}

//====================================================================================================================
//
// important functions
//
//====================================================================================================================
double HeI_Dnem_for_i_effective(double z, int ires){ return calc_HeI_Dnem_i_splines_effective(z, ires); }
double HeI_pd_for_i_effective(double z, int ires){ return calc_HeI_pd_i_splines_effective(z, ires); }

//====================================================================================================================
//
// different Integrals over the spectral distortion
//
//====================================================================================================================
#include "./Solve_PDEs_HeI_integrals.cpp"

//====================================================================================================================
// different Grids
//====================================================================================================================
#include "./Solve_PDEs_HeI_grid.cpp"

//====================================================================================================================
// for xmgrace output
//====================================================================================================================
#ifdef GRACE_DEFINED
#include "./Solve_PDEs_HeI_Grace.cpp"
#endif

//====================================================================================================================
//
// functions to setup PDESolver; this should only happen the first time
//
//====================================================================================================================
PDE_solver_functions_HeI PDE_funcs_HeI;
PDE_Stepper_Data PDE_D_HeI;
bool PDE_funcs_are_set_HeI=0;

//====================================================================================================================
void arm_PDE_solver_HeI(double nu21, Cosmos &cos, Gas_of_HeI_Atoms &HeIAtoms, 
                        vector<int> &lines, vector<double> &xarr, vector<double> &resonances)
{
    if(!PDE_funcs_are_set_HeI)
    {
        //==================================================================
        // here time is wasted a bit, but this is only done once...
        //==================================================================
        PDE_funcs_HeI.HI_Xe=calc_HI_Xe;
        PDE_funcs_HeI.HI_rho=calc_HI_rho;
        PDE_funcs_HeI.HI_X1s=calc_HI_X1s;
        //
        PDE_funcs_HeI.HeI_X1s=calc_HeI_X1s;
        PDE_funcs_HeI.HeI_Xi=calc_HeI_Xi;
        PDE_funcs_HeI.HeI_Dnem_i=HeI_Dnem_for_i_effective;
        PDE_funcs_HeI.HeI_pd_i=HeI_pd_for_i_effective;
        
        Setup_PDE_Data_HeI(cos, HeIAtoms, PDE_funcs_HeI);
        
        PDE_funcs_HeI.res_index_of_lines=lines;
        
        //==================================================================
        int npts=xarr.size();
        
        init_profile_memory_HeI(npts, resonances.size(), PDE_funcs_HeI);
        
        //==================================================================
        // setup for main run
        //==================================================================
        init_PDE_Stepper_Data(PDE_D_HeI, npts);
        if(solve_with_Kernel>0) setup_Lagrange_interpolation_coefficients_O2_Int(PDE_D_HeI, xarr);
        else setup_Lagrange_interpolation_coefficients_O2(PDE_D_HeI, xarr);
        
        //==================================================================
        // HI Lyc photon-ionization cross section at grid-points
        //==================================================================
        for(int k=0; k<(int)xarr.size(); k++)
            PDE_funcs_HeI.sig_Lyc_x[k]=HeIAtoms.HILyc.sig_phot_ion(xarr[k]*nu21);
        
        PDE_funcs_HeI.nu2s1s=HeIAtoms.Sing.Level(2, 0).Get_Dnu_1s2();
        
        unsigned long jem=0;
        locate_JC(&xarr[0], xarr.size(), PDE_funcs_HeI.nu2s1s, &jem);
        PDE_funcs_HeI.emission_2s1s=jem+1;
        
        
        PDE_funcs_are_set_HeI=1;
    }
    
    reset_PDE_solver_variables();
    
    return;
}

//====================================================================================================================
struct Solve_PDE_Data_HeI
{
    int npts;
    int iresmax;
    vector<double> resonances;
    vector<double> Fwork;
    vector<double> xarr, yarr;
};

Solve_PDE_Data_HeI SPDE_D_HeI;
bool SPDE_D_is_set_HeI=0;

//====================================================================================================================
void init_Solve_PDE_Data_HeI(int iresmax, vector<int> &lines, Gas_of_HeI_Atoms &HeIAtoms)
{
    //==================================================================
    if(SPDE_D_is_set_HeI && iresmax!=SPDE_D_HeI.iresmax)
    { 
        cout << " init_Solve_PDE_Data :: the number of resonances has changed."
             << " Please check your code... Resetting. " << endl; 
        
        SPDE_D_is_set_HeI=0;
        PDE_funcs_are_set_HeI=0;
    } 
    
    //==================================================================
    double nu21=HeIAtoms.nP_S_profile(2).Get_nu21();
    double nui1;
    
    if(!SPDE_D_is_set_HeI)
    {
        SPDE_D_HeI.resonances.clear();
        SPDE_D_HeI.Fwork.clear();
        SPDE_D_HeI.xarr.clear(); SPDE_D_HeI.yarr.clear();
        
        //==============================================================
        // set resonances
        //==============================================================
        SPDE_D_HeI.iresmax=iresmax;
        
        for(int n=0; n<=iresmax; n++)
        {
            if(lines[n]==-10)
            { 
                // if the upper point is not used as boundary 
                // condition, just enlarge frequency range.
                SPDE_D_HeI.resonances.push_back(SPDE_D_HeI.resonances.back()*1.05);
                break;
            }
                
            int i_HeI=get_HeI_index(lines[n]);
            
            nui1=HeIAtoms.Get_nu21(HeIAtoms.Get_n(i_HeI), HeIAtoms.Get_l(i_HeI), 
                                   HeIAtoms.Get_S(i_HeI), HeIAtoms.Get_J(i_HeI), 
                                   1, 0, 0, 0);
             
            SPDE_D_HeI.resonances.push_back(nui1/nu21);
        }
        
        //==============================================================
        // prepare grid
        //============================================================== 
        double xmin=xmin_glob_HeI; 
        double xmax=SPDE_D_HeI.resonances.back();
        if(lines[iresmax]==-10) SPDE_D_HeI.resonances.pop_back();

        init_PDE_xarr_cores_HeI(SPDE_D_HeI.xarr, xmin, xmax, SPDE_D_HeI.resonances);

        //==============================================================
        // allocate remaining memory
        //============================================================== 
        SPDE_D_HeI.npts=SPDE_D_HeI.xarr.size();
        SPDE_D_HeI.Fwork.resize(SPDE_D_HeI.npts);
        SPDE_D_HeI.yarr.resize(SPDE_D_HeI.npts);
        
        if(show_messages_HeI>=2) 
            cout << " # of grid-points per resonance: " 
                 << (SPDE_D_HeI.npts-npts_2s1s_HeI)/SPDE_D_HeI.resonances.size()
                 << " Total: " << SPDE_D_HeI.npts << endl;
        
        SPDE_D_is_set_HeI=1;
    }   

    //==================================================================
    // initial solution
    //==================================================================
    for(int k=0; k<(int)SPDE_D_HeI.yarr.size(); k++) SPDE_D_HeI.yarr[k]=0.0;

    return;
}

//====================================================================================================================
//
// to output solution for spectrum
//
//====================================================================================================================
void output_solution_HeI(int output_count, string endname, double nu21, 
                         Cosmos &cos, Gas_of_HeI_Atoms &HeIAtoms, 
                         double zout, int iresmax, bool boundary_up)
{
    string fname=HeI_outputdir+"Sol/sol.HeI."+int_to_string(output_count, 5)+endname;
    ofstream ofile(fname.c_str());
    ofile.precision(8);
    
    //==============================================================================================
    // reference spectrum (at zout)
    //==============================================================================================
    double x_c=const_h_kb*nu21/cos.TCMB(zout);
    double rho=calc_HI_rho(zout);
    double Tm=rho*cos.TCMB(zout);
    //
    double exp_xc=exp(x_c);
    double exp_xc_rho=exp(x_c/rho);
    //
    for(int i=0; i<SPDE_D_HeI.npts; i++)
    { 
        double exp_fac_Tg=PDE_funcs_HeI.exp_x[i]*exp_xc;
        double exp_fac_Te=PDE_funcs_HeI.exp_x[i]*exp_xc_rho;
        
        ofile << SPDE_D_HeI.xarr[i] << " " 
              << HeIAtoms.nP_S_profile(2).nu2x(SPDE_D_HeI.xarr[i]*nu21, Tm) << " " 
              << SPDE_D_HeI.yarr[i] << " " 
              << SPDE_D_HeI.yarr[i]*pow(SPDE_D_HeI.xarr[i], 3) << " ";
                
        ofile << SPDE_D_HeI.yarr[i] << " " 
              << SPDE_D_HeI.yarr[i]/exp_fac_Tg << " " 
              << SPDE_D_HeI.yarr[i]/exp_fac_Te << " "; 
        
        ofile << zout << endl;
    }           
    
    ofile.close();  
    
    return;
}

//====================================================================================================================
void output_distortion_HeI(string fname, double nu21, double z)
{
    ofstream ofile(fname.c_str());
    ofile.precision(8);
    
    ofile << "# HeI spectral distortion at z= " << z 
    << "\n# Columns are:"
    << "\n# x=nu/nu21, y=Dn(x)=nu21*Dn(nu), x^3 y, nu(z), DI(z), nu(z=0) in GHz, DI(z=0)" 
    << "\n# (nu in GHz and DI in ergs cm^-2 sr^-1 Hz^-1 sec^-1)"
    << "\n#" << endl;
    
    double fac=2.0*const_h*const_cl*pow(nu21/const_cl, 3);
    
    for(int i=0; i<SPDE_D_HeI.npts; i++)
    { 
        ofile << SPDE_D_HeI.xarr[i] << " " << SPDE_D_HeI.yarr[i] << " " 
              << SPDE_D_HeI.yarr[i]*pow(SPDE_D_HeI.xarr[i], 3) << " " 
              // spectral distortion at zout
              << SPDE_D_HeI.xarr[i]*nu21*1.0e-9 << " "
              << SPDE_D_HeI.yarr[i]/nu21*fac*pow(SPDE_D_HeI.xarr[i], 3) << " "
              // spectral distortion at z=0
              << SPDE_D_HeI.xarr[i]*nu21*1.0e-9/(1.0+z) << " " 
              << SPDE_D_HeI.yarr[i]/nu21*fac*pow(SPDE_D_HeI.xarr[i]/(1.0+z), 3) 
              << endl;
    }           
    
    ofile.close();  
    
    return;
}


//====================================================================================================================
void output_distortion_HeI(int output_count, string endname, double nu21, double z)
{
    string fname=HeI_outputdir+"Sol/sol.HeI_distortion."+int_to_string(output_count, 5)+endname;
    output_distortion_HeI(fname, nu21, z);
    
    return;
}

//====================================================================================================================
//
// code with effective rates
//
//====================================================================================================================
int compute_HeI_DPesc_with_diffusion_equation_effective(int nHeI_Trans, vector<double> &DF_vec_z, 
                                                        vector<vector<double> > &DF_vec,
                                                        vector<double> &DI1_2s_vec, 
                                                        vector<int> &HeI_level_indices,
                                                        double zs, double ze, 
                                                        Cosmos &cos, 
                                                        Gas_of_HeI_Atoms &HeIAtoms,
                                                        const vector<vector<double> > &HeI_Solution,
                                                        int it_num, bool last_it)
{
    if(nHeI_Trans<2){ cerr << " Not enough helium shells. Exiting. " << endl; exit(0); }

    Set_Load_Populations_HeI_verbosity(show_messages_HeI);
    Set_HeI_pd_Rp_splines_effective_verbosity(show_messages_HeI);    

    //==============================================================================================
    // extension for filename
    //==============================================================================================
    string exten="";
    
    if(solve_with_Kernel==1) exten+=".Kernel_Int";

    //==============================================================================================
	// output for Grace (optional)
    //==============================================================================================
#ifdef GRACE_DEFINED
	bool out_put_GRACE=0;
    bool write_to_disk_GRACE=0;
	int GR_output_step=1;
#endif

    //==============================================================================================
    // general flags
    //==============================================================================================
    bool do_output=0;
    bool write_solution_z0=0;
    bool output_DP=0;
    bool output_DI1=0;
    
    //==============================================================================================
    bool boundary_up=0;                          // change upper boundary condition (==0 --> Dn_x=0)
    
    //==============================================================================================
    // step size
    //==============================================================================================
    double dz_out=10.0;
    int output_count=0;
    int output_step=10;
    int DP_step=1;
    
    // the stability of the iteration depends quite a bit on this, so careful!
    //if(it_num==0){ dz_out=10.0; DP_step=1; }
    //if(it_num==1){ dz_out=5.0; DP_step=2; }
    //if(it_num==2){ dz_out=2.5; DP_step=4; }

    //==============================================================================================
    // frequency rescaling x==nu/nu21
    //==============================================================================================
    double nu21=HeIAtoms.nP_S_profile(2).Get_nu21();
    
    //==============================================================================================
    // select lines
    //==============================================================================================
    HeI_level_indices.clear();

    if(nHeI_Trans>=2)
    {
        HeI_level_indices.push_back(HeIAtoms.Get_Level_index(2, 1, 1, 1)); // 2^3 P - 1 S
        HeI_level_indices.push_back(HeIAtoms.Get_Level_index(2, 1, 0, 1)); // 2^1 P - 1 S
    }
    
    if(nHeI_Trans>=3)
    {
        HeI_level_indices.push_back(HeIAtoms.Get_Level_index(3, 1, 1, 1)); // 3^3 P - 1 S
        HeI_level_indices.push_back(HeIAtoms.Get_Level_index(3, 2, 0, 1)); // 3^1 D - 1 S
        HeI_level_indices.push_back(HeIAtoms.Get_Level_index(3, 1, 0, 1)); // 3^1 P - 1 S
    }
    
    if(nHeI_Trans>=4)
    {
        int nnm=(int) min(10, nHeI_Trans);
        
        for(int nn=4; nn<=nnm; nn++)
        {
            HeI_level_indices.push_back(HeIAtoms.Get_Level_index(nn, 1, 1, 1)); // n^3 P - 1 S
            HeI_level_indices.push_back(HeIAtoms.Get_Level_index(nn, 2, 0, 1)); // n^1 D - 1 S
            HeI_level_indices.push_back(HeIAtoms.Get_Level_index(nn, 1, 0, 1)); // n^1 P - 1 S
        }
    }    
    //==============================================================================================
//    HeI_level_indices.clear();
//    HeI_level_indices.push_back(HeIAtoms.Get_Level_index(3, 1, 0, 1)); 
//    HeI_level_indices.push_back(HeIAtoms.Get_Level_index(2, 1, 1, 1)); 

//    HeI_level_indices.push_back(HeIAtoms.Get_Level_index(2, 1, 1, 1)); // 2^3 P - 1 S
//    HeI_level_indices.push_back(HeIAtoms.Get_Level_index(2, 1, 0, 1)); // 2^1 P - 1 S
    
    
    //==============================================================================================
    // save resolved level indices for computation
    //==============================================================================================
    vector<int> lines;
    for(int line=0; line<(int)HeI_level_indices.size(); line++) 
        lines.push_back(get_res_index_HeI(HeI_level_indices[line]));
    
    if(boundary_up==0) lines.push_back(-10);

    //==============================================================================================
    // total number of treated resonances; the last resonance is 
    // only boundary condition or otherwise is not used.
    //==============================================================================================
    int iresmax=lines.size()-1; 
    
    //==============================================================================================
    // prepare memory for back-communication
    //==============================================================================================
    DF_vec_z.clear();
    DF_vec.clear();
    DI1_2s_vec.clear();
    
    vector<double> DPesc_dum_z;
    for(int line=0; line<(int)HeI_level_indices.size(); line++) DF_vec.push_back(DPesc_dum_z); 
    
    //==============================================================================================
    // load precomputed solutions (load could be avoided...)
    //==============================================================================================
    compute_Xi_HeI_splines(zs, ze, HeIAtoms, HeI_Solution);
    
    set_up_splines_for_HeI_pd_Rp_effective(ze, zs, cos, HeIAtoms, HeI_Solution);
    
    //==============================================================================================
    // memory for Solve_PDE
    //==============================================================================================
    init_Solve_PDE_Data_HeI(iresmax, lines, HeIAtoms);
    
    //==============================================================================================
    // memory for PDE Solver
    //==============================================================================================
    arm_PDE_solver_HeI(nu21, cos, HeIAtoms, lines, SPDE_D_HeI.xarr, SPDE_D_HeI.resonances);
            
    //==============================================================================================
    // outputs/names/etc
    //==============================================================================================
    string fname, endname=exten+".it_"+int_to_string(it_num)+".dat";
    ofstream DPfile, Ifile;
    
    if(output_DP)
    {
        fname=HeI_outputdir+"DPesc/DP.HeI"+endname;
        DPfile.open(fname.c_str());
        DPfile.precision(8);
    }

    if(output_DI1 && HeI_2s_1s_corr_is_on)
    {
        fname=HeI_outputdir+"DPesc/DI1.HeI.2s1s"+endname;
        Ifile.open(fname.c_str());
        Ifile.precision(8);
 
    }
    
    //==============================================================================================
    // main run
    //==============================================================================================
    if(show_messages_HeI>=1) cout << "\n entering HeI-diffusion part " << endl << endl;

    //==============================================================================================
	// prepare grace output
    //==============================================================================================
#ifdef GRACE_DEFINED
	if(out_put_GRACE) prepare_Grace_HeI(0.5, SPDE_D_HeI.resonances, write_to_disk_GRACE);	
#endif
    
    double zin=zs, zout, dz, theta_HeI_loc=theta_HeI;
    int solve_with_Kernel_loc=last_it ? solve_with_Kernel : 0;
    int switched_off_kernel_treatment=100;
    do
    {
        //=======================================================================
        // theta-parameter (added March 2011)
        // Begining with a theta~1 run stabilizes the solution at the beginning,
        // however, highest accuracy in the time-step is reached for theta ~ 0.5
        // which should be used after the initial 'burn in' phase.
        //======================================================================= 
        if(zin>zs-100.0) theta_HeI_loc=0.999;
        else theta_HeI_loc=theta_HeI;
        //theta_HeI_loc=theta_HeI;

        //=======================================================================
        // always add some small amount of diffusion to stabilize the PDE
        //=======================================================================
        set_electron_scattering_efficiency_HeI(1.0);   // full electron scattering with FP-approach
        //set_electron_scattering_efficiency_HeI(0.002); // 0.01 --> 10% difference for Trip;
                                                       // 0.001 --> very small difference; 
        
        //===================================================================
        // switch kernel treatment off at very low Xe_HeI
        //=================================================================== 
        if(solve_with_Kernel_loc==1)
        {
            if(fabs(calc_HeI_X1s(zin)/cos.fHe()-1.0)<=1.0e-3)
            { 
                solve_with_Kernel_loc=0; 
                if(show_messages_HeI>=2) 
                    cout << " switching kernel treatment off ( " << zin << " )" << endl;
                
                switched_off_kernel_treatment=0;
            }
        }
        
        //==========================================================================================
        // set step-size
        //==========================================================================================
        if(switched_off_kernel_treatment<=30) dz=-dz_out/10.0;
        else dz=-dz_out;
        
        zout=max(ze, zin+dz);
        switched_off_kernel_treatment++;

        //==========================================================================================
        // define lower & upper boundary
        //==========================================================================================
        double xeval=SPDE_D_HeI.xarr[0]*(1.0+zin)/(1.0+zout); // (zin > zout)
        int i;
        for(i=0; i<SPDE_D_HeI.npts; i++) if(SPDE_D_HeI.xarr[i]>=xeval){ i--; break; }
        double y_lower, dydum;
        polint_JC(&SPDE_D_HeI.xarr[0], &SPDE_D_HeI.yarr[0], SPDE_D_HeI.npts, xeval, i, 6, &y_lower, &dydum);
        double y_upper=( boundary_up ? PDE_funcs_HeI.HeI_Dnem_i(zout, lines[iresmax])*nu21 : 0.0);
        
        //==========================================================================================
        // solve PDE
        //==========================================================================================
        if(solve_with_Kernel_loc==0) Step_PDE_O2t(theta_HeI_loc, zin, zout, SPDE_D_HeI.xarr, SPDE_D_HeI.yarr,
                                                  y_lower, y_upper, PDE_D_HeI, def_PDE_HeI);       
        
        else if(solve_with_Kernel_loc==1) Step_PDE_O2t_Int(theta_HeI_loc, 1, zin, zout, SPDE_D_HeI.xarr, 
                                                           SPDE_D_HeI.yarr, y_lower, y_upper, PDE_D_HeI, 
                                                           def_PDE_HeI, electron_kernel_corr_func);   
                    
        //==========================================================================================
        // some messages if wanted
        //==========================================================================================
        if(show_messages_HeI>=3) 
            cout << " " << zs << " --> " << zout << " # " << output_count 
                 << " y_low= " << y_lower << " i= " << i << endl;

        //==========================================================================================
        // write results only after some initial evolution
        //==========================================================================================
        if(zout<zs-20.0)
        {
            //======================================================================================
			// output to grace (JPEG)
            //======================================================================================
#ifdef GRACE_DEFINED
            if(solve_with_Kernel_loc>1) 
            {
                compute_kernel_factors(3200, nu21, cos.Te(zout), SPDE_D_HeI.xarr,
                                       PDE_funcs_HeI.Kernel_plus, PDE_funcs_HeI.Kernel_minus, 1.0e-5);
            
                electron_kernel_corr_func(zout, SPDE_D_HeI.xarr, SPDE_D_HeI.yarr, 
                                          PDE_funcs_HeI.Kernel_minus);

                if (!(output_count%GR_output_step) && out_put_GRACE && GraceIsOpen()) 
                    output_Grace_HeI(zout, write_to_disk_GRACE, SPDE_D_HeI.xarr, 
                                     SPDE_D_HeI.yarr, PDE_funcs_HeI.Kernel_plus, 
                                     PDE_funcs_HeI.Kernel_minus, exten, output_count);
            }
            else 
            {
                if (!(output_count%GR_output_step) && out_put_GRACE && GraceIsOpen()) 
                    output_Grace_HeI(zout, write_to_disk_GRACE, SPDE_D_HeI.xarr, 
                                     SPDE_D_HeI.yarr, exten, output_count);
            }
#endif
            //======================================================================================
            // output solution
            //======================================================================================
            if(!(output_count%output_step) && do_output)
            {
                output_solution_HeI(output_count, endname, nu21, cos, HeIAtoms, 
                                    zout, iresmax, boundary_up);
                
                output_distortion_HeI(output_count, endname, nu21, zout);
            }
 
            //======================================================================================
            // compute correction
            //======================================================================================
            if(!(output_count%DP_step))
            {
                //==================================================================================
                DF_vec_z.push_back(zout);
                
                //==================================================================================
                // HeI singlet 2s-1s channel correction
                //==================================================================================
                if(HeI_2s_1s_corr_is_on)    
                    compute_DI1_2s_and_dump_it_HeI(zout, nu21, SPDE_D_HeI.xarr, SPDE_D_HeI.yarr, 
                                                   HeIAtoms, cos, Ifile, DI1_2s_vec, 
                                                   SPDE_D_HeI.Fwork, PDE_funcs_HeI.exp_x,
                                                   phi_2s1s_HeI, PDE_funcs_HeI, output_DI1);
                
                //==================================================================================
                // SPDE_D_HeI.yarr is Dn_x
                //==================================================================================
                for(int line=0; line<(int)HeI_level_indices.size(); line++)
                {
                    
                    DF_vec[line].push_back(compute_DP_HeI(line, zout, nu21, cos, HeIAtoms, 
                                                          SPDE_D_HeI.xarr, SPDE_D_HeI.yarr, 
                                                          SPDE_D_HeI.Fwork, 
                                                          PDE_funcs_HeI.exp_x, PDE_funcs_HeI));
                }

                if(output_DP)
                {
                    double fac_z=calc_HeI_X1s(zout)*cos.NH(zout)/8.0/PI/cos.H(zout);
                    
                    DPfile << zout << " ";
                    for(int line=0; line<(int)DF_vec.size(); line++)
                    {
                        double w=HeIAtoms.Get_gw(HeI_level_indices[line])/HeIAtoms.Get_gw(0);
                        double eta_S=w*phi_i_Ai1s_HeI(line)*pow(const_cl/phi_i_nui1s_HeI(line), 3); 
                        double PS=p_ij(eta_S*fac_z);

                        DPfile << DF_vec[line].back() << " " 
                               << DF_vec[line].back()+PS << " " 
                               << PS << " ";
                    }
                    
                    DPfile << endl;
                }
            }
            
            output_count++;
            
            if(show_messages_HeI>=2) cout << endl;
        }
        
        zin=zout;
    }
    while(zout>ze);
   
    //==============================================================================================
    // write solution for spectrum at z=0
    //==============================================================================================
    if(write_solution_z0)
    {
        fname=HeI_outputdir+"Sol/sol.HeI.z_0"+endname;
        output_distortion_HeI(fname, nu21, ze);
    }
        
    if(output_DI1 && HeI_2s_1s_corr_is_on) Ifile.close();
    if(output_DP) DPfile.close();
    
    //==============================================================================================
	// Flush the output buffer and close Grace 
    //==============================================================================================
#ifdef GRACE_DEFINED
	if(out_put_GRACE && GraceIsOpen()) GraceClose();
#endif
        
    return 0;
}

//====================================================================================================================
//====================================================================================================================
