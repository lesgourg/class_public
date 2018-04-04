//==================================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: compute the first few Raman-profiles
// last modification: Aug 2012
//==================================================================================================
// 18th Aug, 2012: Mnr is now saved to disc after first setup. After that it is only loaded.

#ifndef RAMAN_PROFILES_H
#define RAMAN_PROFILES_H

#include <string>

namespace Raman_profiles 
{
    //==================================================
    // set verbosity level
    //==================================================
    void set_verbosity_Raman(int v);

    void test_Raman_stuff();    
    void init_Raman_profiles();

    //==================================================
    // ns-1s & nd-1s Raman profile functions
    //==================================================
    double sigma_ns_1s_Raman(int n, double y);
    double sigma_nd_1s_Raman(int n, double y);
    double sigma_ns_1s_Raman_ratio(int n, double y);
    double sigma_nd_1s_Raman_ratio(int n, double y);
    void dump_ns_1s_Raman_profile(int ni, string fname);
    void dump_nd_1s_Raman_profile(int ni, string fname);
    
    //==================================================
    // 2s-1s Raman profile functions
    //==================================================
    double sigma_2s_1s_non_res(double y);
    double sigma_2s_1s_res(double y);
    double sigma_2s_1s_res(double y, int choice);   
    double sigma_2s_1s_poles(double y);
    //
    double sigma_2s_1s_Raman(double y);
    double sigma_2s_1s_Raman_ratio(double y);
    //
    void dump_2s_1s_Raman_profile(string fname);

    //==================================================
    // 3s-1s Raman profile functions
    //==================================================
    double sigma_3s_1s_non_res(double y);
    double sigma_3s_1s_res(double y);
    double sigma_3s_1s_res(double y, int choice);   
    double sigma_3s_1s_poles(double y);
    //
    double sigma_3s_1s_Raman(double y);
    double sigma_3s_1s_Raman_ratio(double y);
    //
    void dump_3s_1s_Raman_profile(string fname);

    //==================================================
    // 3d-1s Raman profile functions
    //==================================================
    double sigma_3d_1s_non_res(double y);
    double sigma_3d_1s_res(double y);
    double sigma_3d_1s_res(double y, int choice);   
    double sigma_3d_1s_poles(double y);
    //
    double sigma_3d_1s_Raman(double y);
    double sigma_3d_1s_Raman_ratio(double y);
    //
    void dump_3d_1s_Raman_profile(string fname);
}

#endif
