//==================================================================================================
// Author Jens Chluba
//
// first implementation: July 2010
// last modification   : July 2014
//==================================================================================================
// 30.07.2014: Xi_Data is passed on read only
// 04.06.2011: added nD-1s transitions to transfer problem
// 26.03.2011: added Ly-n integrals & electron scattering Kernel support

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <sstream>

//==================================================================================================
// libs and routines
//==================================================================================================
#include "physical_consts.h"
#include "routines.h"
#include "Cosmos.h"

#include "Sobolev.h"
#include "Patterson.h"

#include "Load.Populations.HI.h"
#include "HI_pd_Rp_splines_effective.h"

#include "PDE_solver.h"
#include "define_PDE.h"
#include "Solve_PDEs.h"

#include "HI_Transition_Data.h"
#include "Raman_profiles.h"
#include "nsnd_2gamma_profiles.h"


//==================================================================================================
using namespace std;
using namespace HI_Transition_Data;
using namespace Raman_profiles;
using namespace nsnd_2gamma_profiles;

string HI_outputdir=COSMORECDIR+"./temp/";

//==================================================================================================
//
// flags for message output
//
//==================================================================================================
int show_messages=-1;


//==================================================================================================
//
// settings for the grid
//
//==================================================================================================
int core_HI=25;            // points per Doppler width in the core region; 
                           // For +/- 10 xD core=25 implies 500 points (which is recommended).
int dec_HI=50;             // points per decade; from 10^-4 --> 1 for ndec==50 gives ~ 200 points

double theta_HI=0.55;

int npts_2s1s=400;         // number of points for 2s-1s two-photon part (300-400 seems enough)
double xmin_glob=0.0001;   // this boundary seems important. A smaller value is better...
//double xmin_glob=0.005;    // previous setting

//==================================================================================================
//
// switch on/off electron scattering kernel corrections
//
//==================================================================================================
int solve_with_Kernel_HI=0; // ==0: no kernel correction; 
                            // ==1: kernel correction using iteration; This runmode takes very long, 
                            //      and is only included for completeness. The additional correction 
                            //      is tiny. This runmode is not recommended.


//==================================================================================================
//
// setting different paramters for PDE problem regarding the HeI Singlet 2s-1s channel
//
//==================================================================================================
bool compute_full_2s_1s_corr=1;
bool HI_2s_1s_corr_is_on=1;

//==================================================================================================
// switch on/off HI 2s-1s channel; this overwrites all other settings for 2s-1s channel
//==================================================================================================
void switch_off_1s_2s_correction()
{ 
    HI_switch_off_2s1s(); 
    HI_2s_1s_corr_is_on=0; 
    
    return; 
}

void switch_on_1s_2s_correction()
{ 
    HI_switch_on_2s1s(); 
    HI_2s_1s_corr_is_on=1; 
    
    return; 
}

//==================================================================================================
// HI 1s-2s absorption can be switched off separately with this.
//==================================================================================================
void switch_off_1s_2s_absorption_correction()
{
    compute_full_2s_1s_corr=0;
    return;
}


//==================================================================================================
//
// important functions
//
//==================================================================================================
double HI_Dnem_for_nl_effective(double z, int n, int l)
{ return calc_HI_Dnem_nl_splines_effective(z, n, l); }

double HI_pd_for_nl_effective(double z, int n, int l)
{ return calc_HI_pd_nl_splines_effective(z, n, l); }

//==================================================================================================
double Dn_nu_ref_effective(int Lyn, double z, double Dn_ref, double pd, double nu, 
                           Cosmos &cos, Gas_of_Atoms &HIA)
{
    double Tm=calc_HI_rho(z)*cos.TCMB(z);
    double tau_d=pd*3.0*HIA.HI_Lyn_profile(Lyn).Get_A21()
                   *pow(HIA.HI_Lyn_profile(Lyn).Get_lambda21(), 3)
                   *calc_HI_X1s(z)*cos.NH(z)/8.0/PI/cos.H(z);
    
    double xi_I=HIA.HI_Lyn_profile(Lyn).xi_Int(HIA.HI_Lyn_profile(Lyn).nu2x(nu, Tm), 1.0e+5, 
                                               HIA.HI_Lyn_profile(Lyn).aVoigt(Tm));
    return Dn_ref*(1.0-exp(-tau_d*xi_I));
}

double Dn_nu_ref_effective_sum(int Lyn_max, double z, double nu, double nu21, 
                               Cosmos &cos, Gas_of_Atoms &HIA)
{
    double r=0.0;
    
    for(int ni=2; ni<=Lyn_max; ni++) 
        r+=Dn_nu_ref_effective(ni, z, HI_Dnem_for_nl_effective(z, ni, 1)*nu21, 
                               HI_pd_for_nl_effective(z, ni, 1), nu, cos, HIA);
    
    return r;
}



//==================================================================================================
//
// different Integrals over the spectral distortion
//
//==================================================================================================
#include "./Solve_PDEs_integrals.cpp"

//==================================================================================================
// different Grids
//==================================================================================================
#include "./Solve_PDEs_grid.cpp"

//==================================================================================================
// for xmgrace output
//==================================================================================================
#ifdef GRACE_DEFINED
#include "./Solve_PDEs_Grace.cpp"
#endif

//==================================================================================================
//
// functions to setup PDESolver; this should only happen the first time
//
//==================================================================================================
PDE_solver_functions PDE_funcs;
PDE_Stepper_Data PDE_D;
bool PDE_funcs_are_set=0;

//==================================================================================================
void arm_PDE_solver(Cosmos &cos, Gas_of_Atoms &HIA, vector<double> &xarr, 
                    vector<double> &resonances)
{
    if(!PDE_funcs_are_set)
    {
        //=========================================================================== 
        // here time is wasted a bit, but this is only done once...
        //=========================================================================== 
        PDE_funcs.HI_Xe=calc_HI_Xe;
        PDE_funcs.HI_rho=calc_HI_rho;
        PDE_funcs.HI_X1s=calc_HI_X1s;
        //
        PDE_funcs.HI_Xnl=calc_HI_Xnl;
        PDE_funcs.HI_Dnem_nl=HI_Dnem_for_nl_effective;
        PDE_funcs.HI_pd_nl=HI_pd_for_nl_effective;
        
        Setup_PDE_Data(cos, HIA, PDE_funcs);
        
        //=========================================================================== 
        int npts=xarr.size();
        int nresmax=resonances.size()+1;
        
        //=========================================================================== 
        // save ratio of 2g & Raman profiles once grid is set
        //=========================================================================== 
        set_verbosity_Raman(show_messages+1);
        init_Raman_profiles();
        //=========================================================================== 
        //dump_ns_1s_Raman_profile(2, COSMORECDIR+"./temp/test.2s.R.dat");
        //
        //for(int n=3; n<=Get_nmax_Raman_correction(); n++)  
        //{
        //    dump_ns_1s_Raman_profile(n, COSMORECDIR+"./temp/test."+int_to_string(n)+"s.R.dat");
        //    dump_nd_1s_Raman_profile(n, COSMORECDIR+"./temp/test."+int_to_string(n)+"d.R.dat");
        //}
        
        set_verbosity_2gamma(show_messages+1);
        init_nsnd_2gamma_profiles();
        //=========================================================================== 
        //dump_ns_1s_2gamma_profile(2, COSMORECDIR+"./temp/test.2s.2g.dat");
        //
        //for(int n=3; n<=Get_nmax_two_g_correction(); n++)  
        //{
        //    dump_ns_1s_2gamma_profile(n, COSMORECDIR+"./temp/test."+int_to_string(n)+"s.2g.dat");
        //    dump_nd_1s_2gamma_profile(n, COSMORECDIR+"./temp/test."+int_to_string(n)+"d.2g.dat");
        //}
        
        init_profile_memory(npts, nresmax, PDE_funcs);
        
        //=========================================================================== 
        // setting profile ratios
        //=========================================================================== 
        double xres_loc;

        //=========================================================================== 
        // Raman profiles
        //=========================================================================== 
        if(Get_nmax_Raman_correction()>=2) 
            for(int k=0; k<npts; k++) 
                PDE_funcs.Raman_ns_1s[0].ratio[k]=sigma_2s_1s_Raman_ratio(xarr[k]-1.0);

        for(int ni=3; ni<=Get_nmax_Raman_correction(); ni++)
        {
            xres_loc=(1.0-1.0/ni/ni)/0.75;
            
            for(int k=0; k<npts; k++)
            {
                double xs=xarr[k]/xres_loc-1.0;
                PDE_funcs.Raman_ns_1s[ni-2].ratio[k]=sigma_ns_1s_Raman_ratio(ni, xs);
                PDE_funcs.Raman_nd_1s[ni-2].ratio[k]=sigma_nd_1s_Raman_ratio(ni, xs);
            }
        }

        //=========================================================================== 
        for(int ni=(int)max(2, Get_nmax_Raman_correction()+1); ni<=nresmax; ni++)
            for(int k=0; k<npts; k++)
                PDE_funcs.Raman_ns_1s[ni-2].ratio[k]=PDE_funcs.Raman_nd_1s[ni-2].ratio[k]=1.0;
        
        //=========================================================================== 
        // two-photon profiles
        //=========================================================================== 
        if(Get_nmax_two_g_correction()>=2) 
            for(int k=0; k<npts; k++) 
                PDE_funcs.two_g_ns_1s[0].ratio[k]=sigma_2s_1s_2gamma(xarr[k]);
        
        for(int ni=3; ni<=Get_nmax_two_g_correction(); ni++)
        {
            xres_loc=(1.0-1.0/ni/ni)/0.75;
        
            for(int k=0; k<npts; k++)
            {
                double xs=xarr[k]/xres_loc;
                PDE_funcs.two_g_ns_1s[ni-2].ratio[k]=sigma_ns_1s_2gamma_ratio(ni, xs);
                PDE_funcs.two_g_nd_1s[ni-2].ratio[k]=sigma_nd_1s_2gamma_ratio(ni, xs);
            }
        }

        //=========================================================================== 
        for(int ni=(int)max(2, Get_nmax_two_g_correction()+1); ni<=nresmax; ni++)
            for(int k=0; k<npts; k++)
                PDE_funcs.two_g_ns_1s[ni-2].ratio[k]=PDE_funcs.two_g_nd_1s[ni-2].ratio[k]=1.0;
        
        //=========================================================================== 
        // A-coefficients and splitting points
        //=========================================================================== 
        int index=0;
        PDE_funcs.index_emission.clear();
        for(; index<npts; index++) if(xarr[index]>=0.5){ PDE_funcs.index_2=index; break; }
        //
        index=0;
        for(int k=0; k<(int)resonances.size(); k++)
        {
            for(int m=0; m<(int)resonances.size(); m++)
            {
                PDE_funcs.Voigt_profiles_A_npns[m][k]=Get_A_npks(m+2, k+2);
                PDE_funcs.Voigt_profiles_A_npnd[m][k]=Get_A_npkd(m+2, k+2);
            }
            //
            for(; index<npts; index++) if(xarr[index]>=resonances[k])
            { PDE_funcs.index_emission.push_back(index); break; }
        }   
        
        if(PDE_funcs.index_emission.size()==0) PDE_funcs.index_emission.push_back(npts);
        
        //=========================================================================== 
        // setup for main run
        //=========================================================================== 
        init_PDE_Stepper_Data(PDE_D, npts);

        if(solve_with_Kernel_HI==0)
            setup_Lagrange_interpolation_coefficients_O2(PDE_D, xarr);
        else if(solve_with_Kernel_HI==1)
            setup_Lagrange_interpolation_coefficients_O2_Int(PDE_D, xarr);
    
        PDE_funcs_are_set=1;
    }
    
    //===========================================================================
    // make sure normalization of A2s1s is reset (added Sept 4th 2014)
    //===========================================================================
    if(Get_nmax_two_g_correction()>=2 && PDE_funcs_are_set)
        for(int k=0; k<(int)xarr.size(); k++)
            PDE_funcs.two_g_ns_1s[0].ratio[k]=sigma_2s_1s_2gamma(xarr[k]);

    
    reset_PDE_solver_variables();
    
    return;
}

//==================================================================================================
struct Solve_PDE_Data
{
    int npts;
    int nresmax;
    vector<double> resonances;
    vector<double> Fwork;
    vector<double> xarr, yarr;
};

Solve_PDE_Data SPDE_D;
bool SPDE_D_is_set=0;

//==================================================================================================
void init_Solve_PDE_Data(int nresmax, bool bound)
{
    //==================================================================
    if(SPDE_D_is_set && nresmax!=SPDE_D.nresmax)
    { 
        cout << " init_Solve_PDE_Data :: the number of resonances has changed."
             << " Please check your code... Resetting. " << endl; 
        
        SPDE_D_is_set=0;
        PDE_funcs_are_set=0;
    } 
    
    //==================================================================
    if(!SPDE_D_is_set)
    {
        SPDE_D.resonances.clear();
        SPDE_D.Fwork.clear();
        SPDE_D.xarr.clear(); SPDE_D.yarr.clear();
        
        //==============================================================
        // set resonances
        //==============================================================
        SPDE_D.nresmax=nresmax;
        for(int n=2; n<=nresmax; n++) 
            SPDE_D.resonances.push_back((1.0-1.0/n/n)/0.75);
        
        //==============================================================
        // prepare grid
        //============================================================== 
        double xmin=xmin_glob; 
        double xmax=( 1.0-1.0/(nresmax+1)/(nresmax+1) )/0.75;
  
        init_PDE_xarr_cores_HI(SPDE_D.xarr, xmin, xmax, SPDE_D.resonances);
        
        //==============================================================
        // allocate remaining memory
        //============================================================== 
        SPDE_D.npts=SPDE_D.xarr.size();
        SPDE_D.Fwork.resize(SPDE_D.npts);
        SPDE_D.yarr.resize(SPDE_D.npts);

        if(show_messages>=2) 
            cout << " # of grid-points per resonance: " 
                 << (SPDE_D.npts-npts_2s1s)/SPDE_D.resonances.size() 
                 << " Total: " << SPDE_D.npts << endl;
        
        SPDE_D_is_set=1;
    }   
    
    //==================================================================
    // initial solution
    //==================================================================
    for(int k=0; k<(int)SPDE_D.yarr.size(); k++) SPDE_D.yarr[k]=0.0;
    
    return;
}

//==================================================================================================
//
// to output solution for spectrum
//
//==================================================================================================
void output_solution(int output_count, string endname, double nu21, Cosmos &cos, Gas_of_Atoms &HIA, 
                     double zout, int nresmax, bool boundary_up)
{
    string fname=HI_outputdir+"Sol/sol.HI."+int_to_string(output_count, 5)+endname;
    ofstream ofile(fname.c_str());
    ofile.precision(8);
    
    fname=HI_outputdir+"Sol/sol.HI.approx."+int_to_string(output_count, 5)+endname;
    ofstream ofileappr(fname.c_str());
    ofileappr.precision(8);
    
    //=========================================================================== 
    // reference spectrum (at zout)
    //=========================================================================== 
    double x_c=const_h_kb*nu21/cos.TCMB(zout);
    double Dnem=HI_Dnem_for_nl_effective(zout, 2, 1);
    double DnL_ref=Dnem*nu21;
    double rho=calc_HI_rho(zout);
    double Tm=rho*cos.TCMB(zout);
    int nmax=(boundary_up ? nresmax+1 : nresmax);
    //
    double exp_xc=exp(x_c);
    double exp_xc_rho=exp(x_c/rho);
    double exp_xc_2=exp(x_c*32.0/27.0);
    //
    ofile << "# z= " << zout << endl;
    //
    for(int i=0; i<SPDE_D.npts; i++)
    { 
        double exp_fac_Tg=PDE_funcs.exp_x[i]*exp_xc;
        double exp_fac_Tg_2=PDE_funcs.exp_x[i]*exp_xc_2;
        double exp_fac_Te=PDE_funcs.exp_x[i]*exp_xc_rho;
        
        ofile << SPDE_D.xarr[i] << " " 
              << HIA.HI_Lyn_profile(2).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " 
              << HIA.HI_Lyn_profile(3).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " 
              << SPDE_D.yarr[i] << " " << SPDE_D.yarr[i]*pow(SPDE_D.xarr[i], 3) << " " 
              << DnL_ref << " " << DnL_ref*exp_fac_Te << " ";
        
        double Dnu_2=Dn_nu_ref_effective(2, zout, HI_Dnem_for_nl_effective(zout, 2, 1)*nu21, 
                                         HI_pd_for_nl_effective(zout, 2, 1), 
                                         SPDE_D.xarr[i]*nu21, cos, HIA);
        
        double Dnu_3=Dn_nu_ref_effective(3, zout, HI_Dnem_for_nl_effective(zout, 3, 1)*nu21, 
                                         HI_pd_for_nl_effective(zout, 3, 1), 
                                         SPDE_D.xarr[i]*nu21, cos, HIA);
        
        double Dnu_tot=Dn_nu_ref_effective_sum(nmax, zout, SPDE_D.xarr[i]*nu21, nu21, cos, HIA);
        
        ofileappr << SPDE_D.xarr[i] << " " 
                  << HIA.HI_Lyn_profile(2).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " 
                  << HIA.HI_Lyn_profile(3).nu2x(SPDE_D.xarr[i]*nu21, Tm) << " " 
        //
        << Dnu_2 << " " << Dnu_3 << " " << Dnu_tot << " "
        //
        << Dnu_2*pow(SPDE_D.xarr[i], 3) << " " 
        << Dnu_3*pow(SPDE_D.xarr[i], 3) << " " 
        << Dnu_tot*pow(SPDE_D.xarr[i], 3) << " "
        //
        << Dnu_2*exp_fac_Tg << " " << Dnu_3*exp_fac_Tg_2 << endl;
        
        ofile << SPDE_D.yarr[i]/exp_fac_Te - DnL_ref << " " 
              << SPDE_D.yarr[i] << " " << SPDE_D.yarr[i]/exp_fac_Tg << " " 
              << SPDE_D.yarr[i]/exp_fac_Te << " "; 
        
        //=======================================================
        // change in the occupation number
        //=======================================================
        ofile << SPDE_D.yarr[i]/nu21 << " ";
        
        //=======================================================
        // Blackbody + distortion
        //=======================================================
        ofile << SPDE_D.yarr[i]/nu21 + PDE_funcs.exp_x[i] << " ";
        
        ofile << endl;
    }           
    
    ofile.close();  
    ofileappr.close();
    
    return;
}

//==================================================================================================
void output_distortion(string fname, double nu21, double z)
{
    ofstream ofile(fname.c_str());
    ofile.precision(8);
    
    ofile << "# HI spectral distortion at z= " << z 
    << "\n# Columns are:"
    << "\n# x=nu/nu21, y=Dn(x)=nu21*Dn(nu), x^3 y, nu(z), DI(z), nu(z=0) in GHz, DI(z=0)" 
    << "\n# (nu in GHz and DI in ergs cm^-2 sr^-1 Hz^-1 sec^-1)"
    << "\n#" << endl;
    
    double fac=2.0*const_h*const_cl*pow(nu21/const_cl, 3);
    
    for(int i=0; i<SPDE_D.npts; i++)
    { 
        ofile << SPDE_D.xarr[i] << " " << SPDE_D.yarr[i] << " " 
              << SPDE_D.yarr[i]*pow(SPDE_D.xarr[i], 3) << " " 
              // spectral distortion at zout
              << SPDE_D.xarr[i]*nu21*1.0e-9 << " "
              << SPDE_D.yarr[i]/nu21*fac*pow(SPDE_D.xarr[i], 3) << " "
              // spectral distortion at z=0
              << SPDE_D.xarr[i]*nu21*1.0e-9/(1.0+z) << " " 
              << SPDE_D.yarr[i]/nu21*fac*pow(SPDE_D.xarr[i]/(1.0+z), 3) 
              << endl;
    }           
    
    ofile.close();  
    
    return;
}


//==================================================================================================
void output_distortion(int output_count, string endname, double nu21, double z)
{
    string fname=HI_outputdir+"Sol/sol.distortion."+int_to_string(output_count, 5)+endname;
    output_distortion(fname, nu21, z);
    
    return;
}

//==================================================================================================
//
// code with effective rates
//
//==================================================================================================
int compute_DPesc_with_diffusion_equation_effective(vector<double> &DF_vec_z, 
                                                    vector<vector<double> > &DF_Ly_n, 
                                                    vector<vector<double> > &DF_2_gamma_vec, 
                                                    vector<vector<double> > &DF_Raman_vec, 
                                                    vector<vector<double> > &DF_nD1s_vec, 
                                                    vector<double> &DI1_2s_vec, 
                                                    double zs, double ze, 
                                                    int nmax_2g_corrs, int nmax_R_corrs,
                                                    Cosmos &cos, Gas_of_Atoms &HIA, 
                                                    const vector<vector<double> > &HI_Solution,
                                                    int it_num)
{
    Set_Load_Populations_HI_verbosity(show_messages);
    Set_HI_pd_Rp_splines_effective_verbosity(show_messages);
    
    //=========================================================================== 
    // initial settings for PDE-problem
    //=========================================================================== 
    HI_switch_on_line_scattering();
    HI_switch_on_e_scattering();
    HI_switch_on_line_em_abs();
    HI_switch_on_2s1s();
    
    //=========================================================================== 
    // extension for filename
    //=========================================================================== 
    string exten="";

    if(solve_with_Kernel_HI==1) exten+=".Kernel_Int";
    
    //=========================================================================== 
	// output for Grace
	//=========================================================================== 
#ifdef GRACE_DEFINED
	bool out_put_GRACE=1;
    bool write_to_disk_GRACE=0;
	int GR_output_step=1;
#endif
    
    //=========================================================================== 
    // general flags
    //=========================================================================== 
    bool do_output=0;
    bool output_DI1=0;
    bool output_2gR=0;
    bool write_solution_z0=0;

    //=========================================================================== 
    // additional correction terms
    //=========================================================================== 
    bool include_Ly_corr=0;   // ==1: include the correction from the 1+1 photon 
                              // terms in the Ly-series. 
    bool include_nD1s_corr=0; // ==1: include the correction from the 1+1 photon 
                              // terms in the nD-1s quadrupole lines. 

    //=========================================================================== 
    int nShells=HIA.Get_nShells();
    bool boundary_up=0;           // change upper boundary condition
    
    //=========================================================================== 
    // step size
    //=========================================================================== 
    double dz_out=10.0;
    int output_count=0;
    int output_step=10;
    int DF_step=1;
    
    // the stability of the iteration depends quite a bit on this, so careful!
    //if(it_num==0){ dz_out=10.0; DF_step=1; }
    //if(it_num==1){ dz_out=5.0; DF_step=2; }
    //if(it_num==2){ dz_out=2.5; DF_step=4; }
    
    //=========================================================================== 
    // frequency rescaling
    //=========================================================================== 
    double nu21=HIA.Level(2, 1).Get_Dnu_1s();
    
    //=========================================================================== 
    // switching on/off Raman/two-gamma profile corrections
    //=========================================================================== 
    // always make sure that Raman-processes are not switched on 
    // above two-gamma processes!
    //---------------------------------------------------------------------------
    nmax_R_corrs=(nmax_R_corrs<nmax_2g_corrs ? nmax_R_corrs : nmax_2g_corrs);
    //
    if(nmax_2g_corrs<2) switch_off_two_g_corrections();
    else switch_on_two_g_corrections(nmax_2g_corrs, nShells);
    //
    if(nmax_R_corrs<2) switch_off_Raman_corrections();
    else switch_on_Raman_corrections(nmax_R_corrs, nShells);
    //
    int nresmax=(int)min(8, nShells);
    if(nmax_2g_corrs>=3) nresmax=(int)min(nmax_2g_corrs, nresmax); 

    //=========================================================================== 
    // this is the case for only 2s-1s correction
    //=========================================================================== 
    if(nmax_2g_corrs==2)
    {
        nresmax=2; 
        boundary_up=1; 
        include_Ly_corr=1;
        dz_out=min(dz_out, 2.0);
    }

    //=========================================================================== 
    // for more than 4 shells the step size of the PDE stepper 
    // has to be decreased. The precision of the output should 
    // then be like ~<0.01%.
    //=========================================================================== 
    //if(nresmax>4){ dz_out=min(dz_out, 2.0); DF_step=10; }
    
    //=========================================================================== 
    // prepare memory for back-communication
    //=========================================================================== 
    DF_vec_z.clear();
    DI1_2s_vec.clear();
    DF_Ly_n.clear();
    DF_2_gamma_vec.clear();
    DF_Raman_vec.clear();
    DF_nD1s_vec.clear();
    
    //=========================================================================== 
    // normal Ly-n integrals
    //=========================================================================== 
    if(include_Ly_corr) DF_Ly_n.resize(nresmax-1); 
    
    //=========================================================================== 
    // only two-gamma integrals when the profile correction 
    // is switched on
    //=========================================================================== 
    DF_2_gamma_vec.resize(2*(Get_nmax_two_g_correction()-2));
    
    //=========================================================================== 
    // only include the Raman integrals when the profile correction 
    // is switched on
    //=========================================================================== 
    DF_Raman_vec.resize(2*(Get_nmax_Raman_correction()-2)+1);
    if(Get_nmax_Raman_correction()<2) DF_Raman_vec.clear();
    
    //=========================================================================== 
    // nD-1s quadrupole lines
    //=========================================================================== 
    if(include_nD1s_corr) DF_nD1s_vec.resize(Get_nmax_two_g_correction()-2); 
    
    //=========================================================================== 
    // load precomputed solutions
    //=========================================================================== 
    compute_Xi_HI_splines(zs, ze, HIA, HI_Solution);

    set_up_splines_for_HI_pd_Rp_effective(ze, zs, cos, HIA, HI_Solution);
    
    //=========================================================================== 
    // memory for Solve_PDE
    //=========================================================================== 
    init_Solve_PDE_Data(nresmax, boundary_up);
    
    //===========================================================================
    // memory for PDE Solver
    //=========================================================================== 
    arm_PDE_solver(cos, HIA, SPDE_D.xarr, SPDE_D.resonances);
    
    //=========================================================================== 
    // outputs/names/etc
    //=========================================================================== 
    string fname, endname=exten+".it_"+int_to_string(it_num)+".dat";
    ofstream Ly_n_file, two_g_file, Raman_file, nD1s_file, Ifile;
    
    if(output_2gR)
    {
        if(include_Ly_corr)
        {
            fname=HI_outputdir+"DPesc/DF.Ly_n"+endname;
            Ly_n_file.open(fname.c_str());
            Ly_n_file.precision(8);
        }

        fname=HI_outputdir+"DPesc/DF.2_gamma"+endname;
        two_g_file.open(fname.c_str());
        two_g_file.precision(8);
        
        fname=HI_outputdir+"DPesc/DF.Raman"+endname;
        Raman_file.open(fname.c_str());
        Raman_file.precision(8);

        fname=HI_outputdir+"DPesc/DF.nD1s"+endname;
        nD1s_file.open(fname.c_str());
        nD1s_file.precision(8);
    }
    
    if(output_DI1 && HI_2s_1s_corr_is_on)
    {
        fname=HI_outputdir+"DPesc/DI1.2s1s"+endname;
        Ifile.open(fname.c_str());
        Ifile.precision(8);
    }
    
    //=========================================================================== 
    // main run
    //=========================================================================== 
    if(show_messages>=1) cout << "\n entering HI-diffusion part " << endl << endl;

	//=========================================================================== 
	// prepare grace output
	//=========================================================================== 
#ifdef GRACE_DEFINED
	if(out_put_GRACE) prepare_Grace(SPDE_D.xarr[0], SPDE_D.resonances, write_to_disk_GRACE);	
#endif
    
    double zin=zs, zout, dz, theta_HI_loc=theta_HI;
    do
    {
        //=======================================================================
        // theta-parameter (added March 2011)
        // Begining with a theta~1 run stabilizes the solution at the beginning,
        // however, highest accuracy in the time-step is reached for theta ~ 0.5
        // which should be used after the initial 'burn in' phase.
        //=======================================================================
        if(zin>zs-100.0) theta_HI_loc=0.999;
        else theta_HI_loc=theta_HI;
        //theta_HI_loc=theta_HI;
        
        //=======================================================================
        // add some small amount of diffusion to stabilize the PDE
        //=======================================================================
        if(solve_with_Kernel_HI==0) set_electron_scattering_efficiency(1.0);
        if(solve_with_Kernel_HI==1) set_electron_scattering_efficiency(1.0e-2);
        
        //=======================================================================
        // set step-size
        //=======================================================================
        dz=-dz_out;
        zout=max(ze, zin+dz);

        //=======================================================================
        // define lower & upper boundary
        //=======================================================================
        double xeval=SPDE_D.xarr[0]*(1.0+zin)/(1.0+zout); // (zin > zout)
        int i;
        for(i=0; i<SPDE_D.npts; i++) if(SPDE_D.xarr[i]>=xeval){ i--; break; }
        double y_lower, dydum;
        polint_JC(&SPDE_D.xarr[0], &SPDE_D.yarr[0], SPDE_D.npts, xeval, i, 6, &y_lower, &dydum);
        double y_upper=( boundary_up ? PDE_funcs.HI_Dnem_nl(zout, nresmax+1, 1)*nu21 : 0.0);
        
        //=======================================================================
        // solve PDE
        //=======================================================================
        if(solve_with_Kernel_HI==0) Step_PDE_O2t(theta_HI_loc, zin, zout, SPDE_D.xarr, SPDE_D.yarr, 
                                                 y_lower, y_upper, PDE_D, def_PDE_Lyn_and_2s1s);       
        
        else if(solve_with_Kernel_HI==1) Step_PDE_O2t_Int(theta_HI_loc, 1, zin, zout, 
                                                          SPDE_D.xarr, SPDE_D.yarr, 
                                                          y_lower, y_upper, PDE_D, 
                                                          def_PDE_Lyn_and_2s1s, 
                                                          electron_kernel_corr_func_HI);  
        
        if(show_messages>=3) 
            cout << " " << zs << " --> " << zout << " # " << output_count 
                 << " y_low= " << y_lower << " i= " << i << endl;

        //=======================================================================
        // write results only after some initial evolution
        //=======================================================================
        if(zout<zs-10.0)
        {
			//===================================================================
			// output to grace (JPEG)
			//===================================================================
#ifdef GRACE_DEFINED
			if (!(output_count%GR_output_step) && out_put_GRACE && GraceIsOpen()) 
				output_Grace(zout, write_to_disk_GRACE, SPDE_D.xarr, 
                             SPDE_D.yarr, exten, output_count);
#endif
            
            //===================================================================
            // output solution
            //===================================================================
            if(!(output_count%output_step) && do_output)
                output_solution(output_count, endname, nu21, cos, HIA, zout, nresmax, boundary_up);
 
            //===================================================================
            // compute correction
            //===================================================================
            if(!(output_count%DF_step))
            {
                //===============================================================
                DF_vec_z.push_back(zout);
                
                //===============================================================
                // SPDE_D.yarr is Dn_x
                //===============================================================
                if(HI_2s_1s_corr_is_on)
                    compute_DI1_2s_and_dump_it(zout, SPDE_D.xarr, SPDE_D.yarr, HIA, cos, Ifile, 
                                               DI1_2s_vec, SPDE_D.Fwork, PDE_funcs.exp_x, 
                                               phi_ns1s_2g, PDE_funcs, output_DI1);
                
                //===============================================================
                // Ly-n integrals (added 26.03.2011)
                //===============================================================
                double DF;
                int count=0;
                if(include_Ly_corr)
                    for(int ni=2; ni<=nresmax; ni++)
                    {
                        DF=compute_DF_Ly_n(ni, nresmax, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                           SPDE_D.Fwork, PDE_funcs.exp_x, PDE_funcs);
                        DF_Ly_n[count++].push_back(DF);
                    }
                
                //===============================================================
                // 2-gamma integrals
                //===============================================================
                if(Get_nmax_two_g_correction()>=3)
                {
                    count=0;
                    for(int ni=3; ni<=nresmax; ni++)
                    {
                        DF=compute_DF_2gamma(ni, 0, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                             SPDE_D.Fwork, PDE_funcs.exp_x, phi_ns1s_2g, PDE_funcs);
                        DF_2_gamma_vec[count++].push_back(DF);
                        
                        DF=compute_DF_2gamma(ni, 2, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                             SPDE_D.Fwork, PDE_funcs.exp_x, phi_nd1s_2g, PDE_funcs);
                        DF_2_gamma_vec[count++].push_back(DF);
                    }
                }
                
                //================================================================
                // Raman integrals
                //================================================================
                if(Get_nmax_Raman_correction()>=2)
                {
                    count=0;
                    
                    DF=compute_DF_Raman(2, 0, nresmax, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                        SPDE_D.Fwork, PDE_funcs.exp_x, phi_ns1s_Raman, PDE_funcs);
                    DF_Raman_vec[count++].push_back(DF);
                    
                    for(int ni=3; ni<=Get_nmax_Raman_correction(); ni++)
                    {
                        DF=compute_DF_Raman(ni, 0, nresmax, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                            SPDE_D.Fwork, PDE_funcs.exp_x, phi_ns1s_Raman, PDE_funcs);
                        DF_Raman_vec[count++].push_back(DF);
                        
                        DF=compute_DF_Raman(ni, 2, nresmax, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                            SPDE_D.Fwork, PDE_funcs.exp_x, phi_nd1s_Raman, PDE_funcs);
                        DF_Raman_vec[count++].push_back(DF);
                    }
                }
                
                //===============================================================
                // nD-1s integrals (added 04.06.2011)
                //===============================================================
                if(Get_nmax_two_g_correction()>=3 && include_nD1s_corr)
                {
                    count=0;
                    for(int ni=3; ni<=Get_nmax_two_g_correction(); ni++)
                    {
                        DF=compute_DF_nD1s(ni, nresmax, zout, cos, HIA, SPDE_D.xarr, SPDE_D.yarr, 
                                           SPDE_D.Fwork, PDE_funcs.exp_x, PDE_funcs);
                        DF_nD1s_vec[count++].push_back(DF);
                    }
                }
                
                //===============================================================
                // output correction integrals
                //===============================================================
                if(output_2gR)
                {
                    Ly_n_file << zout << " ";
                    for(int m=0; m<(int)DF_Ly_n.size(); m++) 
                        Ly_n_file << DF_Ly_n[m].back() << " ";
                    Ly_n_file << endl;

                    two_g_file << zout << " ";
                    for(int m=0; m<(int)DF_2_gamma_vec.size(); m++) 
                        two_g_file << DF_2_gamma_vec[m].back() << " ";
                    two_g_file << endl;
                    
                    Raman_file << zout << " ";
                    for(int m=0; m<(int)DF_Raman_vec.size(); m++) 
                        Raman_file << DF_Raman_vec[m].back() << " ";
                    Raman_file << endl;

                    nD1s_file << zout << " ";
                    for(int m=0; m<(int)DF_nD1s_vec.size(); m++) 
                        nD1s_file << DF_nD1s_vec[m].back() << " ";
                    nD1s_file << endl;
                }
            }
            
            output_count++;
            
            if(show_messages>=2) cout << endl;
        }
        
        zin=zout;
    }
    while(zout>ze);
    
    //=====================================================================================
    // write solution for spectrum at z=0
    //=====================================================================================
    if(write_solution_z0)
    {
        fname=HI_outputdir+"Sol/sol.z_0"+endname;
        output_distortion(fname, nu21, ze);
    }
    
    //=====================================================================================
    // clean up
    //=====================================================================================
    if(output_DI1 && HI_2s_1s_corr_is_on) Ifile.close();
    if(output_2gR)
    {
        if(include_Ly_corr) Ly_n_file.close();
        two_g_file.close();
        Raman_file.close();
        nD1s_file.close();
    }
    
    //=====================================================================================
	// Flush the output buffer and close Grace 
    //=====================================================================================
#ifdef GRACE_DEFINED
	if(out_put_GRACE && GraceIsOpen()) GraceClose();
#endif
    
    return 0;
}

//==================================================================================================
//==================================================================================================
