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

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>

#include "physical_consts.h"
#include "routines.h"
#include "Patterson.h"

#include "HI_matrix_elements.h"
#include "HI_Transition_Data.h"
#include "nsnd_2gamma_profiles.h"

using namespace std;
using namespace HI_matrix_elements;
using namespace HI_Transition_Data;

//==================================================================================================
// Atomic Data for HI atom
//==================================================================================================
const double nsnd_2gamma_profiles_C_sd   =9.0/1024.0*pow(const_alpha, 6)*const_cl
                                                    *const_Ry_inf_icm/(1.0+const_me_mp);

const double nsnd_2gamma_profiles_nu_1sc =const_EH_inf_Hz/(1.0+const_me_mp);

//==================================================================================================
// variables
//==================================================================================================
const int nsnd_2gamma_profiles_nmax=5000;
const int nsnd_2gamma_profiles_xpts=500;  // delete the data-files when changing this parameter !
const double nsnd_2gamma_profiles_xmin=5.0e-5;

const string add_mmax=int_to_string(nsnd_2gamma_profiles_nmax);
const string path=COSMORECDIR+"./Development/Line_profiles/two-photon-data/";
const string ns_Mnr_filename=path+"DGamma.Mnr_ns."+add_mmax+".dat";
const string nd_Mnr_filename=path+"DGamma.Mnr_nd."+add_mmax+".dat";

struct nsnd_2gamma_profiles_Data
{
    int xnpt;
    int memindex;
    vector<double> kappa_res;
    //
    vector<double> f_re;
    vector<double> f_im;
    vector<double> yr;
};

vector<nsnd_2gamma_profiles_Data> nsnd_2gamma_profiles_Data_ns_1s(9);
vector<nsnd_2gamma_profiles_Data> nsnd_2gamma_profiles_Data_nd_1s(9);

//==================================================================================================
// local functions
//==================================================================================================
namespace nsnd_2gamma_profiles_local 
{
    double nuij(int n, int np)
    { return nsnd_2gamma_profiles_nu_1sc*(pow(1.0*np,-2) - pow(1.0*n,-2)); } 
    
    //======================================================================================
    // resonance frequencies (n<ni)
    //======================================================================================
    double y_res(int ni, int n){ return -(1.0/ni/ni-1.0/n/n)/(1.0-1.0/ni/ni); }
    
    //======================================================================================
    // energy factor
    //======================================================================================
    double fn_nsnd_2gamma(int ni, int n, int nf, double y)
    { 
        return 1.0/((pow(1.0*ni,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y) 
             + 1.0/((pow(1.0*nf,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y);
    }

    double fn_nsnd_2gamma_cont(int ni, double x, int nf, double y)
    { 
        return 1.0/(( x*x + pow(1.0*ni,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y) 
             + 1.0/(( x*x + pow(1.0*nf,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y);
    }
    
    //======================================================================================
    // energy factors for the resonances
    //======================================================================================
    double Lorentzian(double a, double b){ return a/(a*a+b*b); }
    
    double fn_nsnd_2gamma_r(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double nuscale=nuij(ni, nf);
        
        //==================================================================================
        // comment: the '+' sign is due to change of phase when using 
        // f= 1/(yp+y-i*d) + 1/(ym-y-i*d) instead of f= 1/(yp+y-i*d) - 1/(y-ym-i*d) as in 
        // CS 2009 paper
        //==================================================================================
        return Lorentzian(dy1, Get_Gamma_np(n)/nuscale) + Lorentzian(dy2, Get_Gamma_np(n)/nuscale);
    }

    double fn_nsnd_2gamma_i(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double nuscale=nuij(ni, nf);
        
        return Lorentzian(Get_Gamma_np(n)/nuscale, dy1) - Lorentzian(Get_Gamma_np(n)/nuscale, dy2);
    }
        
    //======================================================================================
    // normalization factor
    //======================================================================================
    double G_ns(int n){ return nsnd_2gamma_profiles_C_sd    *pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    double G_nd(int n){ return nsnd_2gamma_profiles_C_sd/2.5*pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    
    //======================================================================================
    //======================================================================================

    //======================================================================================
    // ns1s and nd1s integrals over free states
    //======================================================================================
    double (*Cj_x_ptr__)(int ni, double y);
    
    double dIntC1sCj_dx(double logx, void *p)
    {
        double *y=(double *) p;
        double x=exp(logx);
        int ni=y[1];
        return x*C1s(x)*Cj_x_ptr__(ni, x)*fn_nsnd_2gamma_cont(ni, x, 1, y[0]);
    }

    double IntC1sC_nsd(double y, int ni, double (*Cj_x)(int ni, double y))
    {
        if(y>=1.0 || y<=0.0) return 0.0;
        
        double r=0.0;
        double a=log(1.0e-8), b=log(0.5), epsrel=1.0e-8, epsabs=1.0e-16;
        double pars[]={y, (double)ni};
        void *p=&pars;
        Cj_x_ptr__=Cj_x;
        
        r =Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(1.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(10.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(100.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        //
        a=b; b=log(500.0); epsabs=max(epsabs, r*epsrel/2.0);
        r+=Integrate_using_Patterson_adaptive(a, b, epsrel, epsabs, dIntC1sCj_dx, p);
        
        return r;
    }

    //======================================================================================
    // general Matrix element ns/d-->np
    //======================================================================================
    double Matrix_element_j_np_1s(int n, double y, int ni, double (*Rj_np)(int ni, int n))
    { return R1snp(n)*Rj_np(ni, n)*fn_nsnd_2gamma(ni, n, 1, y); }
    
    double Matrix_element_j_1s_total_non_res(double y, int ni, 
                                             double (*Rj_np)(int ni, int n), 
                                             double (*Cj)(int ni, double y))
    {
        double r=IntC1sC_nsd(y, ni, Cj);

        for(int nint=nsnd_2gamma_profiles_nmax; nint>=ni; nint--) 
            r+=Matrix_element_j_np_1s(nint, y, ni, Rj_np);
        
        return r;
    }
    
    //======================================================================================
    // Matrix element Mn_ns_1s
    //======================================================================================
    double Matrix_element_ns_1s_total_non_res(int ni, double y)
    { return Matrix_element_j_1s_total_non_res(y, ni, Rksnp, Cks); }

    double Matrix_element_nd_1s_total_non_res(int ni, double y)
    { return Matrix_element_j_1s_total_non_res(y, ni, Rkdnp, Ckd); }

    //======================================================================================
    // Matrix element using the spline functions after setup
    //======================================================================================
    double compute_Mnr_ns_1s(int n, double y)
    { 
        if(y>0.5) return compute_Mnr_ns_1s(n, 1.0-y);

        return calc_spline_JC(max(nsnd_2gamma_profiles_xmin, y),
                              nsnd_2gamma_profiles_Data_ns_1s[n].memindex)/y/(1.0-y);
    }
    
    double compute_Mnr_nd_1s(int n, double y)
    { 
        if(y>0.5) return compute_Mnr_nd_1s(n, 1.0-y);
        
        return calc_spline_JC(max(nsnd_2gamma_profiles_xmin, y), 
                              nsnd_2gamma_profiles_Data_nd_1s[n].memindex)/y/(1.0-y); 
    }
    
    //======================================================================================
    //======================================================================================
}

using namespace nsnd_2gamma_profiles_local;

//==================================================================================================
//==================================================================================================
namespace nsnd_2gamma_profiles 
{

    //======================================================================================
    // set verbosity level
    //======================================================================================
    int verbosity_level=1;
    void set_verbosity_2gamma(int v){ verbosity_level=v; return; }
    
    //======================================================================================
    // testing 
    //======================================================================================
    void test_nsnd_2gamma_stuff()
    {
        wait_f_r();
    }

    //======================================================================================
    // setup routines
    //======================================================================================
    void alloc_memory(int ni, nsnd_2gamma_profiles_Data &Data_i)
    {
        Data_i.kappa_res.clear();
        Data_i.kappa_res.resize(ni, 0.0);
        Data_i.yr.clear();
        Data_i.yr.       resize(ni, 0.0);

        Data_i.f_re.resize(ni);
        Data_i.f_im.resize(ni);
        
        return;
    }
    
    void init_splines(const vector<double> &xarr, vector<double> &yarr, 
                      int ni,
                      nsnd_2gamma_profiles_Data &Data_i, 
                      double (*M_nr_i)(int ni, double y))
    {
        double fac=1.0;
        // to ensure normalization of 2s-1s profile == 2
        if(ni==2) fac=sqrt(1.0/8.224873118283677);
        
        if(yarr.size()==0)
            for(int k=0; k<nsnd_2gamma_profiles_xpts; k++) 
                yarr.push_back(xarr[k]*(1.0-xarr[k])*M_nr_i(ni, xarr[k])*fac);
        
        Data_i.xnpt=nsnd_2gamma_profiles_xpts;
        Data_i.memindex=calc_spline_coeffies_JC(Data_i.xnpt, 
                                                &xarr[0], &yarr[0], 
                                                "nsnd_2gamma_profiles");
        
        return;
    }
    
    void init_mem_and_data(int ni, 
                           nsnd_2gamma_profiles_Data &Data_i, 
                           double (*Rj_np)(int ni, int n), 
                           const vector<double> &xarr, vector<double> &yarr, 
                           double (*M_nr_i)(int ni, double y))
    {
        alloc_memory(ni, Data_i);
        
        for(int n=2; n<ni; n++) 
        {
            Data_i.kappa_res[n]=R1snp(n)*Rj_np(ni, n);
            Data_i.yr[n]=y_res(ni, n);
        }
        
        init_splines(xarr, yarr, ni, Data_i, M_nr_i);
        
        return;
    }
    
    void dump_mess(string mess)
    {
        if(verbosity_level>0) cout << " dumping 2gamma profile for " << mess << endl;
        return;
    }
    
    //======================================================================================
    // i/o of Mnr
    //======================================================================================
    int read_Mnr(double nmin, ifstream &f, vector<vector<double> > &ya)
    {
        double dum;
        
        for(int k=0; k<nsnd_2gamma_profiles_xpts; k++)
        {
            f >> dum;
            
            for(int n=nmin; n<=8; n++)
            {
                f >> dum;
                ya[n].push_back(dum);
            }
            if(!f) return 1;
        }
        
        return 0;
    }

    int read_Mnr(ifstream &fns, vector<vector<double> > &yarr_ns, 
                 ifstream &fnd, vector<vector<double> > &yarr_nd)
    {
        if(read_Mnr(2, fns, yarr_ns) || read_Mnr(3, fnd, yarr_nd))
        {
            for(int n=0; n<=8; n++) 
            {
                yarr_ns[n].clear();
                yarr_nd[n].clear();
            }
            
            return 1;
        }

        return 0;
    }
    
    void export_Mnr(double nmin, const vector<double> &xa, 
                    ofstream &f, const vector<vector<double> > &ya)
    {
        f.precision(16);
        for(int k=0; k<nsnd_2gamma_profiles_xpts; k++)
        {
            f << scientific << xa[k] << " ";
            for(int n=nmin; n<=8; n++) f << scientific << ya[n][k] << " ";
            f << endl;
        }

        return;
    }
    
    void export_Mnr(const vector<double> xa, 
                    ofstream &fns, const vector<vector<double> > &yarr_ns, 
                    ofstream &fnd, const vector<vector<double> > &yarr_nd)
    {
        export_Mnr(2, xa, fns, yarr_ns);
        export_Mnr(3, xa, fnd, yarr_nd);
        return;
    }

    //======================================================================================
    // setup non-resonant part of 2gamma-profile
    //======================================================================================
    void init_nsnd_2gamma_profiles()
    {       
        if(verbosity_level>0)
            cout << "\n init_nsnd_2gamma_profiles::" 
                 << " Initializing ns-1s and nd-1s 2gamma-profiles up to nmax= 8." 
                 << "\n (the very first time might take a moment) "
                 << endl;

        //======================================================================
        // setup splines for non-resonant part
        //======================================================================
        vector<double> xarr(nsnd_2gamma_profiles_xpts);
        vector<vector<double> > yarr_ns(9);
        vector<vector<double> > yarr_nd(9);
                
        init_xarr(nsnd_2gamma_profiles_xmin, 0.5, &xarr[0], 
                  nsnd_2gamma_profiles_xpts, 0, 0);

        //======================================================================
        // check if transition data is there
        //======================================================================
        bool dump_non_res_M=1;
        ifstream ifilens(ns_Mnr_filename.c_str());
        ifstream ifilend(nd_Mnr_filename.c_str());
        
        if(ifilens && ifilend)
            dump_non_res_M=read_Mnr(ifilens, yarr_ns, ifilend, yarr_nd);
        
        ifilens.close();
        ifilend.close();
                
        //======================================================================
        // 2s
        //======================================================================
        init_splines(xarr, yarr_ns[2], 2, 
                     nsnd_2gamma_profiles_Data_ns_1s[2], 
                     Matrix_element_ns_1s_total_non_res);
        
        //======================================================================
        // compute all profiles splines and kappa_n
        //======================================================================
        for(int ni=3; ni<=8; ni++)
        {
            init_mem_and_data(ni, nsnd_2gamma_profiles_Data_ns_1s[ni], Rksnp, 
                              xarr, yarr_ns[ni], 
                              Matrix_element_ns_1s_total_non_res);
            
            init_mem_and_data(ni, nsnd_2gamma_profiles_Data_nd_1s[ni], Rkdnp, 
                              xarr, yarr_nd[ni], 
                              Matrix_element_nd_1s_total_non_res);
        }

        //======================================================================
        // export Mnr matrix elements if required
        //======================================================================
        if(dump_non_res_M)
        {
            ofstream ofilens(ns_Mnr_filename.c_str());
            ofstream ofilend(nd_Mnr_filename.c_str());
            
            export_Mnr(xarr, ofilens, yarr_ns, ofilend, yarr_nd);
            
            ofilens.close();
            ofilend.close();
        }
        
        if(verbosity_level>0) cout << " init_nsnd_2gamma_profiles:: done " << endl;

        return;
    }
    
    //======================================================================================
    //
    // access to different profiles
    //
    //======================================================================================

    //======================================================================================
    //
    // different cross sections ns/d divided by Gnl(!!!)
    //
    //======================================================================================
    double sigma_nsd_1s_res(double y, int ni, int choice, double Mnr, 
                            nsnd_2gamma_profiles_Data &Data_nsd_1s)
    { 
        double r=0.0;

        // prepare fn tables
        for(int n=2; n<ni; n++) 
        {
            Data_nsd_1s.f_re[n]=fn_nsnd_2gamma_r(ni, n, 1, y);
            Data_nsd_1s.f_im[n]=fn_nsnd_2gamma_i(ni, n, 1, y);
        }
        
        // resonance part
        if(choice==0 || choice==1)
            for(int n=2; n<ni; n++) 
            {
                r+=pow(Data_nsd_1s.kappa_res[n], 2)*( pow(Data_nsd_1s.f_re[n], 2) 
                                                    + pow(Data_nsd_1s.f_im[n], 2) );
                
                for(int m=2; m<n; m++) 
                    r+=2.0*Data_nsd_1s.kappa_res[n]*Data_nsd_1s.kappa_res[m]
                        *( Data_nsd_1s.f_re[n]*Data_nsd_1s.f_re[m] 
                         + Data_nsd_1s.f_im[n]*Data_nsd_1s.f_im[m] );
            }
        
        // interference part
        if(choice==0 || choice==2)
        {
            double r_int=0.0;
            for(int n=2; n<ni; n++) r_int+=Data_nsd_1s.kappa_res[n]*Data_nsd_1s.f_re[n];
            
            r+=2.0*Mnr*r_int;
        }
        
        return r*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    double sigma_nsd_1s_poles(double y, int ni, nsnd_2gamma_profiles_Data &Data_nsd_1s)
    { 
        double r=0.0;
        
        // prepare fn tables
        for(int n=2; n<ni; n++) 
        {
            Data_nsd_1s.f_re[n]=fn_nsnd_2gamma_r(ni, n, 1, y);
            Data_nsd_1s.f_im[n]=fn_nsnd_2gamma_i(ni, n, 1, y);
        }
        
        // sum of resonances
        for(int n=2; n<ni; n++) 
            r+=pow(Data_nsd_1s.kappa_res[n], 2)*( pow(Data_nsd_1s.f_re[n], 2) 
                                                + pow(Data_nsd_1s.f_im[n], 2) );
        
        return r*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_nsd_1s_nsnd_2gamma(double y, int ni, double Mnr, 
                                    nsnd_2gamma_profiles_Data &Data_nsd_1s)
    { 
        double r=Mnr*Mnr;
        
        // prepare fn tables
        for(int n=2; n<ni; n++) 
        {
            Data_nsd_1s.f_re[n]=fn_nsnd_2gamma_r(ni, n, 1, y);
            Data_nsd_1s.f_im[n]=fn_nsnd_2gamma_i(ni, n, 1, y);
        }
        
        // resonance part
        for(int n=2; n<ni; n++) 
        {
            r+=pow(Data_nsd_1s.kappa_res[n], 2)*( pow(Data_nsd_1s.f_re[n], 2) 
                                                + pow(Data_nsd_1s.f_im[n], 2) );
            
            for(int m=2; m<n; m++) 
                r+=2.0*Data_nsd_1s.kappa_res[n]*Data_nsd_1s.kappa_res[m]
                     *(Data_nsd_1s.f_re[n]*Data_nsd_1s.f_re[m] 
                      +Data_nsd_1s.f_im[n]*Data_nsd_1s.f_im[m] );
        }
        
        // interference part
        double r_int=0.0;
        for(int n=2; n<ni; n++) r_int+=Data_nsd_1s.kappa_res[n]*Data_nsd_1s.f_re[n];
        r+=2.0*Mnr*r_int;
        
        return r*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    //
    // different cross sections ns && nd
    //
    //======================================================================================
    double sigma_ns_1s_sum_of_Lorentzians(int ni, double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<ni; n++) 
            r+=Get_A_npks(ni, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))
            *( 
               Lorentzian(Get_Gamma_np(n)/nuij(ni, 1), 
                          y-nsnd_2gamma_profiles_Data_ns_1s[ni].yr[n])
              
              +Lorentzian(Get_Gamma_np(n)/nuij(ni, 1), 
                          y-(1.0-nsnd_2gamma_profiles_Data_ns_1s[ni].yr[n])) 
              );
        
        return r/PI;
    }
    
    double sigma_nd_1s_sum_of_Lorentzians(int ni, double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=2; n<ni; n++) 
            r+=Get_A_npkd(ni, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))
            *( 
               Lorentzian(Get_Gamma_np(n)/nuij(ni, 1), 
                          y-nsnd_2gamma_profiles_Data_nd_1s[ni].yr[n])
              
              +Lorentzian(Get_Gamma_np(n)/nuij(ni, 1), 
                          y-(1.0-nsnd_2gamma_profiles_Data_nd_1s[ni].yr[n])) 
              );
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_ns_1s_2gamma(int n, double y)
    { 
        if(n==2) return sigma_2s_1s_2gamma(y);
        
        if(y<0.0 || y>1.0 || n<2 || n>8) return 0.0;
        return G_ns(n)*sigma_nsd_1s_nsnd_2gamma(y, n, compute_Mnr_ns_1s(n, y), 
                                                nsnd_2gamma_profiles_Data_ns_1s[n]);
    }
    
    double sigma_nd_1s_2gamma(int n, double y)
    { 
        if(y<0.0 || y>1.0 || n<3 || n>8) return 0.0;
        return G_nd(n)*sigma_nsd_1s_nsnd_2gamma(y, n, compute_Mnr_nd_1s(n, y), 
                                                nsnd_2gamma_profiles_Data_nd_1s[n]);
    }
    
    double sigma_ns_1s_2gamma_ratio(int n, double y)
    {  
        if(n==2) return 1.0;
        return sigma_ns_1s_2gamma(n, y)/sigma_ns_1s_sum_of_Lorentzians(n, y); 
    }
    
    double sigma_nd_1s_2gamma_ratio(int n, double y)
    { return sigma_nd_1s_2gamma(n, y)/sigma_nd_1s_sum_of_Lorentzians(n, y); }
    
    //======================================================================================
    // plot profile for ns && nd
    //======================================================================================
    void dump_ns_1s_2gamma_profile(int n, string fname)
    {
        dump_mess(int_to_string(n)+"s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=80000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, 1.0, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << sigma_ns_1s_2gamma(n, xarr[k]) << " " 
                  << sigma_ns_1s_sum_of_Lorentzians(n, xarr[k]) << " " 
                  << sigma_ns_1s_2gamma_ratio(n, xarr[k]) 
                  << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_nd_1s_2gamma_profile(int n, string fname)
    {
        dump_mess(int_to_string(n)+"d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int npy=80000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, 1.0, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << sigma_nd_1s_2gamma(n, xarr[k]) << " " 
                  << sigma_nd_1s_sum_of_Lorentzians(n, xarr[k]) << " " 
                  << sigma_nd_1s_2gamma_ratio(n, xarr[k]) 
                  << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections of 2s-1s
    //
    //======================================================================================
    double sigma_2s_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_ns_1s(2, y);
        return G_ns(2)*const_HI_A2s_1s*Mnr*Mnr*pow(y*(1.0-y), 3);
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_2s_1s_2gamma(double y)
    { 
        if(y<0.0 || y>1.0) return 0.0;
        return sigma_2s_1s_non_res(y);
    }
    
    //======================================================================================
    // plot profile for 2s
    //======================================================================================
    void dump_2s_1s_2gamma_profile(string fname)
    {
        dump_mess("2s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(8);
        
        int npy=10000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, 1.0, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << sigma_2s_1s_2gamma(xarr[k]) << " " 
                  << sigma_2s_1s_non_res(xarr[k]) 
                  << endl;

        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 3s
    //
    //======================================================================================
    double sigma_3s_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_ns_1s(3, y);
        return G_ns(3)*Mnr*Mnr*pow(y*(1.0-y), 3);
    }
    
    //======================================================================================
    double sigma_3s_1s_res(double y, int choice)
    { 
        return G_ns(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_ns_1s(3, y), 
                                        nsnd_2gamma_profiles_Data_ns_1s[3]); 
    }
    
    //======================================================================================
    double sigma_3s_1s_poles(double y)
    { return G_ns(3)*sigma_nsd_1s_poles(y, 3, nsnd_2gamma_profiles_Data_ns_1s[3]); }
    
    //======================================================================================
    double sigma_3s_1s_sum_of_Lorentzians(double y)
    { return sigma_ns_1s_sum_of_Lorentzians(3, y); }
    
    //======================================================================================
    double sigma_3s_1s_res(double y){ return sigma_3s_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3s_1s_2gamma(double y){ return sigma_ns_1s_2gamma(3, y); }
    double sigma_3s_1s_2gamma_ratio(double y){ return sigma_ns_1s_2gamma_ratio(3, y); }
    
    //======================================================================================
    // plot profile for 3s
    //======================================================================================
    void dump_3s_1s_2gamma_profile(string fname)
    {
        dump_mess("3s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=3;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, 1.0, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << xarr[k] /( nsnd_2gamma_profiles_Data_ns_1s[3].yr[which_res] ) << " " 
                  << sigma_3s_1s_2gamma(xarr[k]) << " " 
                  << sigma_3s_1s_non_res(xarr[k]) << " " 
                  << sigma_3s_1s_res(xarr[k]) << " " 
                  << sigma_3s_1s_res(xarr[k], 1) << " " 
                  << sigma_3s_1s_res(xarr[k], 2) << " " 
                  << sigma_3s_1s_2gamma(xarr[k])/sigma_3s_1s_poles(xarr[k]) << " " 
                  << sigma_3s_1s_poles(xarr[k]) << " "
                  << sigma_3s_1s_2gamma(xarr[k])/sigma_3s_1s_sum_of_Lorentzians(xarr[k]) << " " 
                  << sigma_3s_1s_sum_of_Lorentzians(xarr[k]) 
                  << endl;
        
        ofile.close();
        
        return;
    }

    //======================================================================================
    //
    // different cross sections 3d
    //
    //======================================================================================
    double sigma_3d_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_nd_1s(3, y);
        return G_nd(3)*Mnr*Mnr*pow(y*(1.0-y), 3);
    }
    
    //======================================================================================
    double sigma_3d_1s_res(double y, int choice)
    { 
        return G_nd(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_nd_1s(3, y), 
                                        nsnd_2gamma_profiles_Data_nd_1s[3]); 
    }
    
    //======================================================================================
    double sigma_3d_1s_poles(double y)
    { return G_nd(3)*sigma_nsd_1s_poles(y, 3, nsnd_2gamma_profiles_Data_nd_1s[3]); }
    
    //======================================================================================
    double sigma_3d_1s_sum_of_Lorentzians(double y)
    { return sigma_nd_1s_sum_of_Lorentzians(3, y); }
    
    //======================================================================================
    double sigma_3d_1s_res(double y){ return sigma_3d_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3d_1s_2gamma(double y){ return sigma_nd_1s_2gamma(3, y); }
    double sigma_3d_1s_2gamma_ratio(double y){ return sigma_nd_1s_2gamma_ratio(3, y); }
    
    //======================================================================================
    // plot profile for 3d
    //======================================================================================
    void dump_3d_1s_2gamma_profile(string fname)
    {
        dump_mess("3d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=3;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, 1.0, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << xarr[k] /( nsnd_2gamma_profiles_Data_nd_1s[3].yr[which_res] ) << " " 
                  << sigma_3d_1s_2gamma(xarr[k]) << " " 
                  << sigma_3d_1s_non_res(xarr[k]) << " " 
                  << sigma_3d_1s_res(xarr[k]) << " " 
                  << sigma_3d_1s_res(xarr[k], 1) << " " 
                  << sigma_3d_1s_res(xarr[k], 2) << " " 
                  << sigma_3d_1s_2gamma(xarr[k])/sigma_3d_1s_poles(xarr[k]) << " " 
                  << sigma_3d_1s_poles(xarr[k]) << " "
                  << sigma_3d_1s_2gamma(xarr[k])/sigma_3d_1s_sum_of_Lorentzians(xarr[k]) << " " 
                  << sigma_3d_1s_sum_of_Lorentzians(xarr[k]) 
                  << endl;
        
        ofile.close();
        
        return;
    }
}

//==================================================================================================
//==================================================================================================
