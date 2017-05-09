//==================================================================================================
// Author Jens Chluba Sept/Oct 2010
// purpose: compute the first few two-photon-profiles
// last modification: July 2014
//==================================================================================================
// 23rd Jul, 2014: 2s-1s profile changed to have normalization set by external const_HI_A2s_1s
// 18th Aug, 2012: Mnr is now saved to disc after first setup. After that it is only loaded.
// 17th Aug, 2012: fixed bug for Rnd->np matrix element + 3d1s profile.
//                 changed spline-setup for non-resonant part to xmin - 0.5;
//                 symmetry is used for x>0.5

#ifndef NSND_2GAMMA_PROFILES_H
#define NSND_2GAMMA_PROFILES_H

#include <string>

namespace nsnd_2gamma_profiles 
{
    //==================================================
    // set verbosity level
    //==================================================
    void set_verbosity_2gamma(int v);
    
    void test_nsnd_2gamma_stuff();  
    void init_nsnd_2gamma_profiles();

    //==================================================
    // ns-1s & nd-1s two-photon profile functions
    //==================================================
    double sigma_ns_1s_2gamma(int n, double y);
    double sigma_nd_1s_2gamma(int n, double y);
    double sigma_ns_1s_2gamma_ratio(int n, double y);
    double sigma_nd_1s_2gamma_ratio(int n, double y);
    void dump_ns_1s_2gamma_profile(int n, string fname);
    void dump_nd_1s_2gamma_profile(int n, string fname);

    //==================================================
    // 2s-1s two-photon profile functions
    //==================================================
    double sigma_2s_1s_non_res(double y);
    double sigma_2s_1s_2gamma(double y);
    void dump_2s_1s_2gamma_profile(string fname);
    
    //==================================================
    // 3s-1s two-photon profile functions
    //==================================================
    double sigma_3s_1s_non_res(double y);
    double sigma_3s_1s_res(double y);
    double sigma_3s_1s_res(double y, int choice);   
    double sigma_3s_1s_poles(double y);
    //
    double sigma_3s_1s_2gamma(double y);
    double sigma_3s_1s_2gamma_ratio(double y);
    void dump_3s_1s_2gamma_profile(string fname);

    //==================================================
    // 3d-1s two-photon profile functions
    //==================================================
    double sigma_3d_1s_non_res(double y);
    double sigma_3d_1s_res(double y);
    double sigma_3d_1s_res(double y, int choice);   
    double sigma_3d_1s_poles(double y);
    //
    double sigma_3d_1s_2gamma(double y);
    double sigma_3d_1s_2gamma_ratio(double y);
    void dump_3d_1s_2gamma_profile(string fname);
}

#endif

//==================================================================================================
//==================================================================================================
