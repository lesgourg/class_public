//==================================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: compute the first few Raman-profiles
// last modification: Aug 2012
//==================================================================================================
// 18th Aug, 2012: Mnr is now saved to disc after first setup. After that it is only loaded.

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
#include "Raman_profiles.h"

using namespace std;
using namespace HI_matrix_elements;
using namespace HI_Transition_Data;

//==================================================================================================
// Atomic Data for HI atom
//==================================================================================================
const double Raman_profiles_C_sd   =9.0/1024.0*pow(const_alpha, 6)*const_cl
                                              *const_Ry_inf_icm/(1.0+const_me_mp);

const double Raman_profiles_nu_1sc =const_EH_inf_Hz/(1.0+const_me_mp);

//==================================================================================================
// variables
//==================================================================================================
const int Raman_profiles_nmax=5000;
const int Raman_profiles_xpts=500;  // delete the data-files when changing this parameter !
const double Raman_profiles_xmin=5.0e-5;

const int Raman_profiles_nres_up=10;  // upper limit on resonances. It will not be possible to  
                                      // get the correct Raman-profile 'above' the last included
                                      // resonance.

const string add_mmax=int_to_string(Raman_profiles_nmax);
const string path=COSMORECDIR+"./Development/Line_profiles/two-photon-data/";
const string ns_Mnr_filename=path+"Raman.Mnr_ns."+add_mmax+".dat";
const string nd_Mnr_filename=path+"Raman.Mnr_nd."+add_mmax+".dat";

struct Raman_profiles_Data
{
    int res_up;
    int xnpt;
    int memindex;
    double xmax;
    vector<double> kappa_res;
    //
    vector<double> f_re;
    vector<double> f_im;
    vector<double> yr;
};

vector<Raman_profiles_Data> Raman_profiles_Data_ns_1s(9);
vector<Raman_profiles_Data> Raman_profiles_Data_nd_1s(9);

//==================================================================================================
// local functions
//==================================================================================================
namespace Raman_profiles_local 
{
    double nuij(int n, int np){ return Raman_profiles_nu_1sc*(pow(1.0*np,-2) - pow(1.0*n,-2)); } 
    
    //======================================================================================
    // resonance frequencies
    //======================================================================================
    double y_res(int ni, int n){ return (1.0/ni/ni-1.0/n/n)/(1.0-1.0/ni/ni); }
    
    //======================================================================================
    // energy factor
    //======================================================================================
    double fn_Raman(int ni, int n, int nf, double y)
    { 
        return 1.0/((pow(1.0*ni,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y) 
             + 1.0/((pow(1.0*nf,-2) -pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y);
    }

    double fn_Raman_cont(int ni, double x, int nf, double y)
    { 
        return 1.0/(( x*x + pow(1.0*ni,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y) 
             + 1.0/(( x*x + pow(1.0*nf,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y);
    }
    
    //======================================================================================
    // energy factors for the resonances
    //======================================================================================
    double Lorentzian(double a, double b){ return a/(a*a+b*b); }
    
    double fn_Raman_r(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double nuscale=nuij(ni, nf);
        
        //==================================================================================
        // comment: the '+' sign is due to change of phase when using 
        // f= 1/(yp+y-i*d) + 1/(ym-y-i*d) instead of f= 1/(yp+y-i*d) - 1/(y-ym-i*d) as in 
        // CS 2009 paper
        //==================================================================================
        return Lorentzian(dy1, Get_Gamma_np(n)/nuscale) + Lorentzian(dy2, Get_Gamma_np(n)/nuscale);
    }

    double fn_Raman_i(int ni, int n, int nf, double y)
    {
        double dy1=(pow(1.0*ni,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) - y;
        double dy2=(pow(1.0*nf,-2)-pow(1.0*n,-2))/(pow(1.0*nf,-2) - pow(1.0*ni,-2)) + y;
        double nuscale=nuij(ni, nf);
        
        return Lorentzian(Get_Gamma_np(n)/nuscale, dy1) - Lorentzian(Get_Gamma_np(n)/nuscale, dy2);
    }
        
    //======================================================================================
    // normalization factor
    //======================================================================================
    double G_ns(int n){ return Raman_profiles_C_sd    *pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    double G_nd(int n){ return Raman_profiles_C_sd/2.5*pow(4.0/3.0*(1.0-1.0/n/n), 5); }
    
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
        return x*C1s(x)*Cj_x_ptr__(ni, x)*fn_Raman_cont(ni, x, 1, y[0]);
    }
    
    
    double IntC1sC_nsd(double y, int ni, double (*Cj_x)(int ni, double))
    {
        // y=nu/nu_n1 --> pole at y=yc
        if(y>=1.0/(ni*ni-1.0)) return 0.0;
        
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
    { return R1snp(n)*Rj_np(ni, n)*fn_Raman(ni, n, 1, y); }
    
    double Matrix_element_j_1s_total_non_res(int nmin, double y, int ni, 
                                             double (*Rj_np)(int ni, int n), 
                                             double (*Cj)(int ni, double y))
    {
        double r=IntC1sC_nsd(y, ni, Cj);
        
        nmin=(int)max(ni+1, nmin);
        
        for(int nint=Raman_profiles_nmax; nint>nmin; nint--) 
            r+=Matrix_element_j_np_1s(nint, y, ni, Rj_np);
        
        for(int nint=ni; nint>=2; nint--)                  
            r+=Matrix_element_j_np_1s(nint, y, ni, Rj_np);

        return r;
    }
    
    //======================================================================================
    // Matrix element Mn_ns_1s
    //======================================================================================
    double Matrix_element_ns_1s_total_non_res(int nmin, int ni, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, ni, Rksnp, Cks); }
    
    double Matrix_element_nd_1s_total_non_res(int nmin, int ni, double y)
    { return Matrix_element_j_1s_total_non_res(nmin, y, ni, Rkdnp, Ckd); }
    
    //======================================================================================
    // Matrix element using the spline functions after setup
    //======================================================================================
    double compute_Mnr_ns_1s(int n, double y)
    { 
        return calc_spline_JC(max(Raman_profiles_xmin, y), 
                              Raman_profiles_Data_ns_1s[n].memindex)/y/(1.0/(n*n-1.0)-y); 
    }
    
    double compute_Mnr_nd_1s(int n, double y)
    { 
        return calc_spline_JC(max(Raman_profiles_xmin, y), 
                              Raman_profiles_Data_nd_1s[n].memindex)/y/(1.0/(n*n-1.0)-y); 
    }

    //======================================================================================
    //======================================================================================
}

using namespace Raman_profiles_local;

//==================================================================================================
//==================================================================================================
namespace Raman_profiles 
{

    //======================================================================================
    // set verbosity level
    //======================================================================================
    int verbosity_level=1;
    void set_verbosity_Raman(int v){ verbosity_level=v; return; }
    
    //======================================================================================
    // testing 
    //======================================================================================
    void test_Raman_stuff()
    {
        cout << C1s(0.000001) << " " << C2s(0.000001) << endl;
        cout << C1s(10.0) << " " << C2s(30) << endl;
        cout << R1snp(10) << " " << R2snp(130) << endl;

        cout << endl;
        cout << R2snp(10) << " " << R2snp(130) << endl;
        cout << R3snp(10) << " " << R3snp(130) << endl;
        cout << R4snp(10) << " " << R4snp(130) << endl;
        cout << R5snp(10) << " " << R5snp(130) << endl;

        cout << endl;
        cout << R3dnp(10) << " " << R3dnp(130) << endl;
        cout << R4dnp(10) << " " << R4dnp(130) << endl;
        cout << R5dnp(10) << " " << R5dnp(130) << endl;
        
        cout << endl;
        cout << C2s(0.00001) << " " << C2s(30) << endl;
        cout << C3s(0.00001) << " " << C3s(30) << endl;
        cout << C4s(0.00001) << " " << C4s(30) << endl;
        cout << C5s(0.00001) << " " << C5s(30) << endl;

        cout << endl;
        cout << C3d(0.00001) << " " << C3d(30) << endl;
        cout << C4d(0.00001) << " " << C4d(30) << endl;
        cout << C5d(0.00001) << " " << C5d(30) << endl;
        
        wait_f_r();
    }

    //======================================================================================
    // setup routines
    //======================================================================================
    void alloc_memory(int nres_up, Raman_profiles_Data &Data_i)
    {
        Data_i.res_up=nres_up;
        Data_i.kappa_res.clear();
        Data_i.kappa_res.resize(Data_i.res_up+1, 0.0);
        Data_i.yr.clear();
        Data_i.yr.       resize(Data_i.res_up+1, 0.0);

        Data_i.f_re.resize(Data_i.res_up+1);
        Data_i.f_im.resize(Data_i.res_up+1);
        
        return;
    }
    
    void init_splines(int ni,
                      const vector<double> &xarr, vector<double> &yarr, 
                      Raman_profiles_Data &Data_i, 
                      double (*M_nr_i)(int nmin, int ni, double y))
    {
        if(yarr.size()==0)
            for(int k=0; k<Raman_profiles_xpts; k++) 
                yarr.push_back(xarr[k]*(1.0/(ni*ni-1.0)-xarr[k])*M_nr_i(Data_i.res_up,ni,xarr[k]));
        
        Data_i.xnpt=Raman_profiles_xpts;
        Data_i.xmax=xarr[Raman_profiles_xpts-1];
        Data_i.memindex=calc_spline_coeffies_JC(Data_i.xnpt, 
                                                &xarr[0], &yarr[0],
                                                "Raman_profiles");
        
        return;
    }
    
    void init_mem_and_data(int ni, int nres_up, 
                           Raman_profiles_Data &Data_i, 
                           double (*Rj_np)(int ni, int n),
                           const vector<double> &xarr, vector<double> &yarr, 
                           double (*M_nr_i)(int nmin, int ni, double y))
    {
        alloc_memory(nres_up, Data_i);
        
        for(int n=ni+1; n<=nres_up; n++) 
        {
            Data_i.kappa_res[n]=R1snp(n)*Rj_np(ni, n);
            Data_i.yr[n]=y_res(ni, n);
        }
        
        //==================================================================================
        // rescale frequency grid for computation of splines
        //==================================================================================
        double x_scale=nuij(nres_up, 1)/nuij(ni, 1)-1.0-Raman_profiles_xmin;
        vector<double> xc(Raman_profiles_xpts);
        for(int k=0; k<Raman_profiles_xpts; k++) xc[k]=xarr[k]*x_scale;
        
        init_splines(ni, xc, yarr, Data_i, M_nr_i);
        
        return;
    }
    
    void dump_mess(string mess)
    {
        if(verbosity_level>0) cout << " dumping Raman profile for " << mess << endl;
        return;
    }
    
    //======================================================================================
    // i/o of Mnr
    //======================================================================================
    int read_Mnr(double nmin, ifstream &f, vector<vector<double> > &ya)
    {
        double dum;
        
        for(int k=0; k<Raman_profiles_xpts; k++)
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
        for(int k=0; k<Raman_profiles_xpts; k++)
        {
            f << scientific << xa[k] << " ";
            for(int n=nmin; n<=8; n++) f << scientific << ya[n][k] << " ";
            f << endl;
        }
        
        return;
    }
    
    void export_Mnr(const vector<double> &xa, 
                    ofstream &fns, const vector<vector<double> > &yarr_ns, 
                    ofstream &fnd, const vector<vector<double> > &yarr_nd)
    {
        export_Mnr(2, xa, fns, yarr_ns);
        export_Mnr(3, xa, fnd, yarr_nd);
        return;
    }

    //======================================================================================
    // setup non-resonant part of Raman-profile; this depends on chosen xmax
    //======================================================================================
    void init_Raman_profiles()
    {       
        if(verbosity_level>0) 
            cout << "\n init_Raman_profiles::" 
                 << " Initializing ns-1s and nd-1s Raman-profiles up to nmax=8 and nres=10." 
                 << "\n (the very first time might take a moment) "
                 << endl;

        //======================================================================
        // setup splines for non-resonant part
        //======================================================================
        vector<double> xarr(Raman_profiles_xpts);
        vector<vector<double> > yarr_ns(9);
        vector<vector<double> > yarr_nd(9);
        
        // all x are scaled to the smale range
        init_xarr(Raman_profiles_xmin, 1.0, &xarr[0], 
                  Raman_profiles_xpts, 0, 0);
    
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
        // ns && nd
        //======================================================================
        init_mem_and_data(2, Raman_profiles_nres_up, 
                          Raman_profiles_Data_ns_1s[2], Rksnp, 
                          xarr, yarr_ns[2], 
                          Matrix_element_ns_1s_total_non_res);
        
        //======================================================================
        // compute all profiles splines and kappa_n
        //======================================================================
        for(int ni=3; ni<=8; ni++)
        {
            init_mem_and_data(ni, Raman_profiles_nres_up, 
                              Raman_profiles_Data_ns_1s[ni], Rksnp, 
                              xarr, yarr_ns[ni], 
                              Matrix_element_ns_1s_total_non_res);

            init_mem_and_data(ni, Raman_profiles_nres_up, 
                              Raman_profiles_Data_nd_1s[ni], Rkdnp, 
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

        if(verbosity_level>0) cout << " init_Raman_profiles:: done " << endl;

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
                            Raman_profiles_Data &Data_nsd_1s)
    { 
        double r=0.0;
        int nmax=(int)min(10, Data_nsd_1s.res_up);
        
        // prepare fn tables
        for(int n=ni+1; n<=nmax; n++) 
        {
            Data_nsd_1s.f_re[n]=fn_Raman_r(ni, n, 1, y);
            Data_nsd_1s.f_im[n]=fn_Raman_i(ni, n, 1, y);
        }
        
        // resonance part
        if(choice==0 || choice==1)
            for(int n=ni+1; n<=nmax; n++) 
            {
                r+=pow(Data_nsd_1s.kappa_res[n], 2)*( pow(Data_nsd_1s.f_re[n], 2) 
                                                    + pow(Data_nsd_1s.f_im[n], 2) );
                
                for(int m=ni+1; m<n; m++) 
                    r+=2.0*Data_nsd_1s.kappa_res[n]*Data_nsd_1s.kappa_res[m]
                        *( Data_nsd_1s.f_re[n]*Data_nsd_1s.f_re[m] 
                         + Data_nsd_1s.f_im[n]*Data_nsd_1s.f_im[m] );
            }
        
        // interference part
        if(choice==0 || choice==2)
        {
            double r_int=0.0;
            for(int n=ni+1; n<=nmax; n++) 
                r_int+=Data_nsd_1s.kappa_res[n]*Data_nsd_1s.f_re[n];
            
            r+=2.0*Mnr*r_int;
        }
        
        return r*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    double sigma_nsd_1s_poles(double y, int ni, Raman_profiles_Data &Data_nsd_1s)
    { 
        double r=0.0;
        int nmax=(int)min(10, Data_nsd_1s.res_up);
        
        // prepare fn tables
        for(int n=ni+1; n<=nmax; n++) 
        {
            Data_nsd_1s.f_re[n]=fn_Raman_r(ni, n, 1, y);
            Data_nsd_1s.f_im[n]=fn_Raman_i(ni, n, 1, y);
        }
        
        // sum of resonances
        for(int n=ni+1; n<=nmax; n++) 
            r+=pow(Data_nsd_1s.kappa_res[n], 2)*( pow(Data_nsd_1s.f_re[n], 2) 
                                                + pow(Data_nsd_1s.f_im[n], 2) );
        
        return r*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_nsd_1s_Raman(double y, int ni, double Mnr, Raman_profiles_Data &Data_nsd_1s)
    { 
        double r=Mnr*Mnr;
        int nmax=(int)min(10, Data_nsd_1s.res_up);
        
        // prepare fn tables
        for(int n=ni+1; n<=nmax; n++) 
        {
            Data_nsd_1s.f_re[n]=fn_Raman_r(ni, n, 1, y);
            Data_nsd_1s.f_im[n]=fn_Raman_i(ni, n, 1, y);
        }
        
        // resonance part
        for(int n=ni+1; n<=nmax; n++) 
        {
            r+=pow(Data_nsd_1s.kappa_res[n], 2)*( pow(Data_nsd_1s.f_re[n], 2) 
                                                + pow(Data_nsd_1s.f_im[n], 2) );
            
            for(int m=ni+1; m<n; m++) 
                r+=2.0*Data_nsd_1s.kappa_res[n]*Data_nsd_1s.kappa_res[m]
                     *(Data_nsd_1s.f_re[n]*Data_nsd_1s.f_re[m] 
                      +Data_nsd_1s.f_im[n]*Data_nsd_1s.f_im[m] );
        }
        
        // interference part
        double r_int=0.0;
        for(int n=ni+1; n<=nmax; n++) r_int+=Data_nsd_1s.kappa_res[n]*Data_nsd_1s.f_re[n];
        r+=2.0*Mnr*r_int;
        
        return r*pow(y*(1.0+y), 3);
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
        for(int n=ni+1; n<=(int)min(10, Raman_profiles_Data_ns_1s[ni].res_up); n++) 
            r+=Get_A_npks(ni, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))
              *Lorentzian(Get_Gamma_np(n)/nuij(ni, 1), y-Raman_profiles_Data_ns_1s[ni].yr[n]);
        
        return r/PI;
    }
    
    double sigma_nd_1s_sum_of_Lorentzians(int ni, double y)
    { 
        double r=0.0;
        
        // sum of resonances
        for(int n=ni+1; n<=(int)min(10, Raman_profiles_Data_nd_1s[ni].res_up); n++) 
            r+=Get_A_npkd(ni, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))
              *Lorentzian(Get_Gamma_np(n)/nuij(ni, 1), y-Raman_profiles_Data_nd_1s[ni].yr[n]);
        
        return r/PI;
    }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_ns_1s_Raman(int ni, double y)
    { 
        if(ni<2 || y<0.0) return 0.0;
        return G_ns(ni)*sigma_nsd_1s_Raman(y, ni, compute_Mnr_ns_1s(ni, y), 
                                           Raman_profiles_Data_ns_1s[ni]);
    }
    double sigma_nd_1s_Raman(int ni, double y)
    { 
        if(ni<3 || y<0.0) return 0.0;
        return G_nd(ni)*sigma_nsd_1s_Raman(y, ni, compute_Mnr_nd_1s(ni, y), 
                                          Raman_profiles_Data_nd_1s[ni]);
    }
    
    double sigma_ns_1s_Raman_ratio(int ni, double y)
    { return sigma_ns_1s_Raman(ni, y)/sigma_ns_1s_sum_of_Lorentzians(ni, y); }
    
    double sigma_nd_1s_Raman_ratio(int ni, double y)
    { return sigma_nd_1s_Raman(ni, y)/sigma_nd_1s_sum_of_Lorentzians(ni, y); }
    
    //======================================================================================
    // plot profile for ns && nd
    //======================================================================================
    void dump_ns_1s_Raman_profile(int ni, string fname)
    {
        dump_mess(int_to_string(ni)+"s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=ni+1;
        int npy=80000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_ns_1s[ni].xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << xarr[k] /( Raman_profiles_Data_ns_1s[ni].yr[which_res] ) << " " 
                  << sigma_ns_1s_Raman(ni, xarr[k]) << " " 
                  << sigma_ns_1s_sum_of_Lorentzians(ni, xarr[k]) << " " 
                  << sigma_ns_1s_Raman_ratio(ni, xarr[k]) 
                  << endl;
        
        ofile.close();
        
        return;
    }
    
    void dump_nd_1s_Raman_profile(int ni, string fname)
    {
        dump_mess(int_to_string(ni)+"d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=ni+1;
        int npy=80000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_nd_1s[ni].xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                 << xarr[k] /( Raman_profiles_Data_nd_1s[ni].yr[which_res] ) << " " 
                 << sigma_nd_1s_Raman(ni, xarr[k]) << " " 
                 << sigma_nd_1s_sum_of_Lorentzians(ni, xarr[k]) << " " 
                 << sigma_nd_1s_Raman_ratio(ni, xarr[k]) 
                 << endl;
        
        ofile.close();
        
        return;
    }
    
    
    //======================================================================================
    //
    // different cross sections 2s
    //
    //======================================================================================
    double sigma_2s_1s_non_res(double y)
    { 
        double Mnr=compute_Mnr_ns_1s(2, y);
        return G_ns(2)*Mnr*Mnr*pow(y*(1.0+y), 3);
    }

    //======================================================================================
    double sigma_2s_1s_res(double y, int choice)
    { 
        return G_ns(2)*sigma_nsd_1s_res(y, 2, choice, compute_Mnr_ns_1s(2, y), 
                                        Raman_profiles_Data_ns_1s[2]); 
    }
    
    //======================================================================================
    double sigma_2s_1s_poles(double y)
    { return G_ns(2)*sigma_nsd_1s_poles(y, 2, Raman_profiles_Data_ns_1s[2]); }

    //======================================================================================
    double sigma_2s_1s_sum_of_Lorentzians(double y)
    { return sigma_ns_1s_sum_of_Lorentzians(2, y); }

    //======================================================================================
    double sigma_2s_1s_res(double y){ return sigma_2s_1s_res(y, 0); }

    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_2s_1s_Raman(double y)
    { 
        if(y<0.0) return 0.0;
        
        return G_ns(2)*sigma_nsd_1s_Raman(y, 2, compute_Mnr_ns_1s(2, y), 
                                          Raman_profiles_Data_ns_1s[2]);
    }
    
    //======================================================================================
    // total cross section with motion
    //======================================================================================
    double sigma_2s_1s_Raman_motion(double y, vector<double> phi_i_y)
    { 
        double phi_V_tot=0.0;
        for(int n=3; n<=(int)min(Raman_profiles_Data_ns_1s[2].res_up, min(10, phi_i_y.size()+3)); n++) 
            phi_V_tot+=Get_A_npks(2, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*phi_i_y[n-3];
            
        return phi_V_tot*sigma_2s_1s_Raman(y)/sigma_2s_1s_sum_of_Lorentzians(y);
    }
    
    double sigma_2s_1s_Raman_ratio(double y)
    { return sigma_ns_1s_Raman_ratio(2, y); }

    //======================================================================================
    // plot profile for 2s
    //======================================================================================
    void dump_2s_1s_Raman_profile(string fname)
    {
        dump_mess("2s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=3;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_ns_1s[2].xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << xarr[k] /( Raman_profiles_Data_ns_1s[2].yr[which_res] ) << " " 
                  << sigma_2s_1s_Raman(xarr[k]) << " " 
                  << sigma_2s_1s_non_res(xarr[k]) << " " 
                  << sigma_2s_1s_res(xarr[k]) << " " 
                  << sigma_2s_1s_res(xarr[k], 1) << " " 
                  << sigma_2s_1s_res(xarr[k], 2) << " " 
                  << sigma_2s_1s_Raman(xarr[k])/sigma_2s_1s_poles(xarr[k]) << " " 
                  << sigma_2s_1s_poles(xarr[k]) << " "
                  << sigma_2s_1s_Raman(xarr[k])/sigma_2s_1s_sum_of_Lorentzians(xarr[k]) << " " 
                  << sigma_2s_1s_sum_of_Lorentzians(xarr[k]) 
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
        return G_ns(3)*Mnr*Mnr*pow(y*(1.0+y), 3);
    }
    
    //======================================================================================
    double sigma_3s_1s_res(double y, int choice)
    { return G_ns(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_ns_1s(3, y), 
                                      Raman_profiles_Data_ns_1s[3]); }
    
    //======================================================================================
    double sigma_3s_1s_poles(double y)
    { return G_ns(3)*sigma_nsd_1s_poles(y, 3, Raman_profiles_Data_ns_1s[3]); }
    
    //======================================================================================
    double sigma_3s_1s_sum_of_Lorentzians(double y)
    { return sigma_ns_1s_sum_of_Lorentzians(3, y); }
    
    //======================================================================================
    double sigma_3s_1s_res(double y){ return sigma_3s_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3s_1s_Raman(double y)
    { return sigma_ns_1s_Raman(3, y); }
    
    //======================================================================================
    // total cross section with motion
    //======================================================================================
    double sigma_3s_1s_Raman_motion(double y, vector<double> phi_i_y)
    { 
        double phi_V_tot=0.0;
        for(int n=4; n<=(int)min(Raman_profiles_Data_ns_1s[3].res_up, min(10, phi_i_y.size()+4)); n++) 
            phi_V_tot+=Get_A_npks(3, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*phi_i_y[n-4];
        
        return phi_V_tot*sigma_3s_1s_Raman(y)/sigma_3s_1s_sum_of_Lorentzians(y);
    }
    
    double sigma_3s_1s_Raman_ratio(double y){ return sigma_ns_1s_Raman_ratio(3, y); }
    
    //======================================================================================
    // plot profile for 3s
    //======================================================================================
    void dump_3s_1s_Raman_profile(string fname)
    {
        dump_mess("3s-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=4;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_ns_1s[3].xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
            << xarr[k] /( Raman_profiles_Data_ns_1s[3].yr[which_res] ) << " " 
            << sigma_3s_1s_Raman(xarr[k]) << " " 
            << sigma_3s_1s_non_res(xarr[k]) << " " 
            << sigma_3s_1s_res(xarr[k]) << " " 
            << sigma_3s_1s_res(xarr[k], 1) << " " 
            << sigma_3s_1s_res(xarr[k], 2) << " " 
            << sigma_3s_1s_Raman(xarr[k])/sigma_3s_1s_poles(xarr[k]) << " " 
            << sigma_3s_1s_poles(xarr[k]) << " "
            << sigma_3s_1s_Raman(xarr[k])/sigma_3s_1s_sum_of_Lorentzians(xarr[k]) << " " 
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
        return G_nd(3)*Mnr*Mnr*pow(y*(1.0+y), 3);
    }
    
    //======================================================================================
    double sigma_3d_1s_res(double y, int choice)
    { 
        return G_nd(3)*sigma_nsd_1s_res(y, 3, choice, compute_Mnr_nd_1s(3, y), 
                                      Raman_profiles_Data_nd_1s[3]); 
    }
    
    //======================================================================================
    double sigma_3d_1s_poles(double y)
    { return G_nd(3)*sigma_nsd_1s_poles(y, 3, Raman_profiles_Data_nd_1s[3]); }
    
    //======================================================================================
    double sigma_3d_1s_sum_of_Lorentzians(double y)
    { return sigma_nd_1s_sum_of_Lorentzians(3, y); }
    
    //======================================================================================
    double sigma_3d_1s_res(double y){ return sigma_3d_1s_res(y, 0); }
    
    //======================================================================================
    // total cross section
    //======================================================================================
    double sigma_3d_1s_Raman(double y)
    { return sigma_nd_1s_Raman(3, y); }
    
    //======================================================================================
    // total cross section with motion
    //======================================================================================
    double sigma_3d_1s_Raman_motion(double y, vector<double> phi_i_y)
    { 
        double phi_V_tot=0.0;
        for(int n=4; n<=(int)min(Raman_profiles_Data_nd_1s[3].res_up, min(10, phi_i_y.size()+4)); n++) 
            phi_V_tot+=Get_A_npkd(3, n)*Get_A_np1s(n)/(FOURPI*Get_Gamma_np(n))*phi_i_y[n-4];
        
        return phi_V_tot*sigma_3d_1s_Raman(y)/sigma_3d_1s_sum_of_Lorentzians(y);
    }
    
    double sigma_3d_1s_Raman_ratio(double y){ return sigma_nd_1s_Raman_ratio(3, y); }
    
    //======================================================================================
    // plot profile for 3d
    //======================================================================================
    void dump_3d_1s_Raman_profile(string fname)
    {
        dump_mess("3d-1s");
        ofstream ofile(fname.c_str());
        ofile.precision(10);
        
        int which_res=4;
        int npy=60000;
        vector<double> xarr(npy);
        init_xarr(1.0e-4, Raman_profiles_Data_nd_1s[3].xmax, &xarr[0], npy, 1, 0);
        
        for(int k=0; k<npy; k++)
            ofile << xarr[k] << " " 
                  << xarr[k] /( Raman_profiles_Data_nd_1s[3].yr[which_res] ) << " " 
                  << sigma_3d_1s_Raman(xarr[k]) << " " 
                  << sigma_3d_1s_non_res(xarr[k]) << " " 
                  << sigma_3d_1s_res(xarr[k]) << " " 
                  << sigma_3d_1s_res(xarr[k], 1) << " " 
                  << sigma_3d_1s_res(xarr[k], 2) << " " 
                  << sigma_3d_1s_Raman(xarr[k])/sigma_3d_1s_poles(xarr[k]) << " " 
                  << sigma_3d_1s_poles(xarr[k]) << " "
                  << sigma_3d_1s_Raman(xarr[k])/sigma_3d_1s_sum_of_Lorentzians(xarr[k]) << " " 
                  << sigma_3d_1s_sum_of_Lorentzians(xarr[k]) 
                  << endl;
        
        ofile.close();
        
        return;
    }
}

//==================================================================================================
//==================================================================================================
