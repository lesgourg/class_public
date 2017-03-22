//==================================================================================================
// accuracy for integrator
//==================================================================================================
double epsrel_HI=1.0e-5;

//#define DIFF_CORR_STOREII   // Changes between two ways of storing the diffusion correction.
                            // Not setting it is consistent with CosmoRec v1.5
                            // make sure to also change this in ./Modules/Diffusion_correction.cpp

//==================================================================================================
// to compute the integrals over the solution for the photon field
//==================================================================================================
struct HI_nx_spline_Data
{
    int memindex;
    double xmin, xmax;
};

bool DF_integrals_splines_are_set=0;
HI_nx_spline_Data Spline_for_nx_integral;

//==================================================================================================
// simple functions for Sobolev approximation
//==================================================================================================
double tau_S_function(int Lyn, double z, double X1s, Cosmos &cos, Gas_of_Atoms &HIA)
{ 
    return 3.0*HIA.HI_Lyn_profile(Lyn).Get_A21()
              *pow(HIA.HI_Lyn_profile(Lyn).Get_lambda21(), 3)
              *X1s*cos.NH(z)/8.0/PI/cos.H(z); 
}

double tau_S_function_nD1s(int n, double z, double X1s, Cosmos &cos, Gas_of_Atoms &HIA)
{ 
    return 5.0*HIA.HI_nD1s_profile(n).Get_A21()
              *pow(HIA.HI_nD1s_profile(n).Get_lambda21(), 3)
              *X1s*cos.NH(z)/8.0/PI/cos.H(z); 
}

//==================================================================================================
//
// functions for integrals of two-gamma & Raman process & remaining 1+1 photon corrections
//
//==================================================================================================
double dnbar_dx_2g_R(double x){ return calc_spline_JC(x, Spline_for_nx_integral.memindex); }

double calc_DF(double xmin, double xmax, double epsabs=1.0e-60)
{
    //=========================================================================
    // Integrals
    //=========================================================================
    double epsrel=epsrel_HI;
    double r=Integrate_using_Patterson_adaptive(xmin, xmax, epsrel, epsabs, dnbar_dx_2g_R);
    
    return r;   
}

double compute_integral_over_resonances(double xmin, double xmax, 
                                        int nmin, int nmax, 
                                        double nu21, Gas_of_Atoms &HIA)
{
    double a=xmin, b;
    double dx_x=1.0e-4;
    double r=0.0, r1, xres;
    double epsabs=1.0e-60;

    //----------------------------------------------------------------------
    // lower boundary may contain large phase space density 
    // (better to separate it)
    //----------------------------------------------------------------------
    b=a*(1.0+0.005); // JC: for nmax=10 this needs to be 0.005 rather than 0.01
    r1=calc_DF(a, b, epsabs);
    r+=r1; epsabs=fabs(r)*epsrel_HI;
    a=b;
    
    for(int n=nmin; n<=nmax; n++)
    {
        xres=HIA.HI_Lyn_profile(n).Get_nu21()/nu21;
        //----------------------------------------------------------------------
        // below Lyn-resonance (term from Sobolev xmin-->0)
        //----------------------------------------------------------------------
        b=xres*(1.0-dx_x);
        r1=calc_DF(a, b, epsabs);
        r+=r1; epsabs=fabs(r)*epsrel_HI;
        a=b;
        
        //----------------------------------------------------------------------
        // across Lyn-resonance
        //----------------------------------------------------------------------
        b=xres*(1.0+dx_x);
        r1=calc_DF(a, b, epsabs);
        r+=r1; epsabs=fabs(r)*epsrel_HI;
        a=b;
    }
    
    b=xmax*(1.0-dx_x);
    r1=calc_DF(a, b, epsabs);
    r+=r1; epsabs=fabs(r)*epsrel_HI;
    
    r1=calc_DF(b, xmax, epsabs);
    r+=r1;

    return r;
} 

//==================================================================================================
//
// remaining 1+1 photon corrections
//
//==================================================================================================
double compute_DF_Ly_n(int ni, int nmax, double z, Cosmos &cos, Gas_of_Atoms &HIA, 
                       vector<double> &xarr, vector<double> &yarr, 
                       vector<double> &Fwork, vector<double> &exp_x, 
                       PDE_solver_functions &PDE_funcs)
{
    if(ni>11){ cerr << " compute_DF_Ly_n:: too many intermediate resonances..." << endl; exit(0); }
    
    int nnu_pts=xarr.size();
    double nui1=HIA.Level(ni, 1).Get_Dnu_1s();
    double nu21=HIA.Level(2 , 1).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    //
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HI_Xnl(z, ni, 1);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s/3.0-exp_i1;
    //
    double Dnem=PDE_funcs.HI_Dnem_nl(z, ni, 1);
    double Dnem_eff=PDE_funcs.Voigt_profiles_Dnem_eff[ni-2];
    double DDnem=Dnem-Dnem_eff;
    double pd=phi_Ly_n_pd(ni);
    double pd_eff=PDE_funcs.Voigt_profiles_pd_eff[ni-2];
    
    //================================================================================
    // Ly-n Sobolev approximation part
    //================================================================================
    double tau_S=tau_S_function(ni, z, X1s, cos, HIA); 
    double PS=p_ij(tau_S), P_d=p_ij(pd*tau_S); 
    double PS_eff=pd*P_d/(1.0-(1.0-pd)*P_d);
    double DRtot=-3.0*X1s*HIA.HI_Lyn_profile(ni).Get_A21()*DnLj*(PS-PS_eff); // no-scattering correction
    DRtot+=pd_eff*3.0*X1s*HIA.HI_Lyn_profile(ni).Get_A21()*(DDnem-Dnem*P_d); // correction to normal Pd
    
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_i1*(1.0/exp_x[k]-1.0); // exp(-xi1)/npl
                
        Fwork[k]=(yarr[k]*f_x-nu21*Dnem_eff)*phi_Ly_n(ni, k); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], 
                                                                &Fwork[0], "DF-integrals");
                                                                
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //==========================================================================
    // compute integral
    //==========================================================================
    double r=compute_integral_over_resonances(xarr[2], xarr[nnu_pts-2], 2, nmax, nu21, HIA);
    
    //==========================================================================
    r*=pd_eff*3.0*X1s*HIA.HI_Lyn_profile(ni).Get_A21()/phi_Ly_n_Dnu(ni); 
    double res=(r-DRtot)/DnLj;
    
    if(show_messages>=2) cout << " difference:: (ni, 1) = " << ni << " Dr-Ly-n = " << r      
                              << " DRtot= " << DRtot << " " << res << endl;
    
    return res;
}

//==================================================================================================
//
// nD-1s 1+1 photon correction (added June 2011)
//
//==================================================================================================
double compute_DF_nD1s(int ni, int nmax, double z, Cosmos &cos, Gas_of_Atoms &HIA, 
                       vector<double> &xarr, vector<double> &yarr, 
                       vector<double> &Fwork, vector<double> &exp_x, 
                       PDE_solver_functions &PDE_funcs)
{
    if(ni>11){ cerr << " compute_DF_Ly_n:: too many intermediate resonances..." << endl; exit(0); }

    Voigtprofile_Dawson *p=&HIA.HI_nD1s_profile(ni);
    
    int nnu_pts=xarr.size();
    double nui1=p->Get_nu21();
    double nu21=HIA.Level(2 , 1).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    double Te=Tg*PDE_funcs.HI_rho(z);
    //
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HI_Xnl(z, ni, 2);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s/5.0-exp_i1;
    double aV  =p->aVoigt(Te);
    
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_i1*(1.0/exp_x[k]-1.0); // exp(-xi1)/npl
        double xD=p->nu2x(nu21*xarr[k], Te);
        
        Fwork[k]=(yarr[k]*f_x-nu21*DnLj)*p->phi(xD, aV); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], 
                                                                &Fwork[0], "DF-integrals");
                                                                
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //==========================================================================
    // compute integral
    //==========================================================================
    double r=compute_integral_over_resonances(xarr[2], xarr[nnu_pts-2], 2, nmax, nu21, HIA);
    
    //==========================================================================
    r*=5.0*X1s*p->Get_A21()/p->DnuT(Te); 
    
    double res=r/DnLj;
    
    if(show_messages>=2) cout << " difference:: (ni, 2) = " << ni << " Dr-nD-1s = " << r      
                              << " " << res << endl;
    
    return res;
}

//==================================================================================================
//
// two-gamma-process
//
//==================================================================================================
double compute_DF_2gamma(int ni, int li, double z, Cosmos &cos, Gas_of_Atoms &HIA, 
                         vector<double> &xarr, vector<double> &yarr, 
                         vector<double> &Fwork, vector<double> &exp_x,
                         double (*phi_2g)(int n, int k), 
                         PDE_solver_functions &PDE_funcs)
{
    if(ni>11){ cerr << " compute_DF_2gamma:: too many intermediate resonances..." << endl; exit(0);}
    
    int nnu_pts=PDE_funcs.index_emission[ni-2];
    double nuik, nui1=HIA.Level(ni, li).Get_Dnu_1s();
    double nu21=HIA.Level(2, 1).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    //
    double w=(2.0*li+1);
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HI_Xnl(z, ni, li);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s/w-exp_i1;
    //
    double Dnem=HI_Dnem_for_nl_effective(z, ni, li);
    double DRtot=0.0;
    
    //================================================================================
    // create data for part from Sobolev approximation
    //================================================================================
    for(int n=2; n<ni; n++) 
    {
        nuik=HIA.Level(ni, li).Get_nu21(n, 1);
        double Aik=HIA.Level(ni, li).Get_A21(n, 1);
        //
        double exp_n=exp(-const_h_kb*HIA.HI_Lyn_profile(n).Get_nu21()/Tg);
        double tau_S_n=tau_S_function(n, z, X1s, cos, HIA);
        //
        double pem=1.0-calc_HI_pd_nl_splines_effective(z, n, 1);
        double exp_ik=exp(-const_h_kb*nuik/Tg);
        double PS=p_ij(tau_S_n);
        double Xkp=calc_HI_Xnl(z, n, 1);
        double Atilde=Aik*pem/(1.0-exp_ik);
        //
        double DnLkp=Xkp/X1s/3.0-exp_n;
        double Del=DnLkp-DnLj/exp_ik;
        //
        DRtot+=Atilde*exp_ik*(Del-PS*DnLkp);
    }
        
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_i1*(1.0/exp_x[k]-1.0); // exp(-xi1)/npl

        Fwork[k]=(yarr[k]*f_x-nu21*Dnem)*phi_2g(ni, k); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], 
                                                                &Fwork[0], "DF-integrals");
                                                                
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //==========================================================================
    // compute integral
    //==========================================================================
    double xmin=max(0.5*nui1/nu21, xarr[0]);
    double xmax=min(nui1/nu21, xarr[nnu_pts-2]);
    double r=compute_integral_over_resonances(xmin, xmax, 2, ni-1, nu21, HIA);
    
    //==========================================================================
    if(show_messages>=2) cout << " difference:: (ni, li)= " << ni << " " << li << " Dr-2g = " << r      
                              << " DRtot= " << DRtot << " " << w*(r-DRtot)/DnLj 
                              << " ( error-check " << DnLj/Dnem-1.0 << " )" << endl;
    
#ifdef DIFF_CORR_STOREII
    return w*X1s*(r-DRtot)/DnLj;
#else
    return w*(r-DRtot);
#endif
}

//==================================================================================================
//
// Raman-process
//
//==================================================================================================
double compute_DF_Raman(int ni, int li, int nmax, double z, Cosmos &cos, Gas_of_Atoms &HIA, 
                        vector<double> &xarr, vector<double> &yarr, 
                        vector<double> &Fwork, vector<double> &exp_x,
                        double (*phi_R)(int n, int k), PDE_solver_functions &PDE_funcs)
{
    if(nmax>10){ cerr<< " compute_DF_Raman:: too many intermediate resonances..." << endl; exit(0);}
    
    int k0=PDE_funcs.index_emission[ni-2];
    int nnu_pts=xarr.size()-k0;
    double nuik, nui1=HIA.Level(ni, li).Get_Dnu_1s();
    double nu21=HIA.Level(2, 1).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    //
    double w=(2.0*li+1.0);
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HI_Xnl(z, ni, li);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s/w-exp_i1;
    //
    double Dnem=HI_Dnem_for_nl_effective(z, ni, li);
    double DRtot=0.0;
    
    //================================================================================
    // create data for part from Sobolev approximation
    //================================================================================
    for(int n=ni+1; n<=nmax; n++) 
    {
        nuik=HIA.Level(n, 1).Get_nu21(ni, li);
        double Aik=HIA.Level(n, 1).Get_A21(ni, li);
        //
        double exp_n=exp(-const_h_kb*HIA.HI_Lyn_profile(n).Get_nu21()/Tg);
        double tau_S_n=tau_S_function(n, z, X1s, cos, HIA);
        //
        double pem=1.0-calc_HI_pd_nl_splines_effective(z, n, 1);
        double exp_ik=exp(-const_h_kb*nuik/Tg);
        double PS=p_ij(tau_S_n);
        double Xkp=calc_HI_Xnl(z, n, 1);
        double Atilde=3.0/w*Aik*pem/(1.0-exp_ik);
        //
        double DnLkp=Xkp/X1s/3.0-exp_n;
        double Del=DnLkp-DnLj*exp_ik;
        //
        DRtot+=Atilde*(Del-PS*DnLkp);
    }
    
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_i1*(1.0/exp_x[k+k0]-1.0); // exp(-xi1)/npl

        Fwork[k+k0]=(yarr[k+k0]*f_x-nu21*Dnem)*phi_R(ni, k+k0); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0+k0], 
                                                                &Fwork[0+k0], "DF-integrals");
                                                                
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, 
                                   &xarr[0+k0], &Fwork[0+k0]);
    
    //==========================================================================
    // compute integral
    //==========================================================================
    double xmin=max(nui1/nu21, xarr[0+k0]);
    double xmax=xarr[nnu_pts+k0-2];
    double r=compute_integral_over_resonances(xmin, xmax, ni+1, nmax, nu21, HIA);
    
    //==========================================================================
    if(show_messages>=2) cout << " difference:: (ni, li)= " << ni << " " << li 
                              << " Dr-Raman = " << r 
                              << " DRtot= " << DRtot << " " << w*(r-DRtot)/DnLj 
                              << " ( error-check " << DnLj/Dnem-1.0 << " )" << endl;
    
#ifdef DIFF_CORR_STOREII
    return w*X1s*(r-DRtot)/DnLj;
#else
    return w*(r-DRtot);
#endif
}


//==================================================================================================
//
// DI1 output function; correction to the 2s-1s channel
//
//==================================================================================================
double calc_DI1_2s_Patt(HI_nx_spline_Data &SD)
{
    //=========================================================================
    // 1s-->2s Integral; only high nu part (nu/nu21>=0.5) matters
    //=========================================================================
    double r=0.0, rel=epsrel_HI, abs=1.0e-60;
    
    r=Integrate_using_Patterson_adaptive(max(SD.xmin, 0.5),min(SD.xmax, 1.0),rel,abs,dnbar_dx_2g_R); 
    
    return r/const_HI_A2s_1s;
}

//==================================================================================================
void compute_DI1_2s_and_dump_it(double z, vector<double> &xarr, vector<double> &yarr, 
                                Gas_of_Atoms &HIA, Cosmos &cos, 
                                ofstream &Ifile, vector<double> &DI1_2s_vec, 
                                vector<double> &Fwork, vector<double> &exp_x,
                                double (*phi_2g)(int n, int k), 
                                PDE_solver_functions &PDE_funcs,
                                bool output=0)
{
    //==========================================================================
    int nnu_pts=PDE_funcs.index_emission[0];
    double nu21=HIA.Level(2, 0).Get_Dnu_1s();
    double Tg=cos.TCMB(z);
    double exp_i1=exp(-const_h_kb*nu21/Tg);
    double Dnem=HI_Dnem_for_nl_effective(z, 2, 0);

    Spline_for_nx_integral.xmin=xarr[0];
    Spline_for_nx_integral.xmax=xarr[nnu_pts-1];
    
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_i1*(1.0/exp_x[k]-1.0); // exp(-xi1)/npl
        
        if(compute_full_2s_1s_corr) Fwork[k]=(yarr[k]*f_x-nu21*Dnem)*phi_2g(2, k); 
        else Fwork[k]=-nu21*Dnem*phi_2g(2, k); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DF_integrals_splines_are_set) 
    {
        Spline_for_nx_integral.memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], 
                                                                &Fwork[0], "DF-integrals");
                                                                
        DF_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral.memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //==========================================================================
    // compute integrals
    //==========================================================================
    double DI1=calc_DI1_2s_Patt(Spline_for_nx_integral);

#ifdef DIFF_CORR_STOREII
    DI1=DI1/Dnem+1.0;
#else
    double Xnl=calc_HI_Xnl(z, 2, 0);
    double X1s=calc_HI_X1s(z);
    double DnLj=Xnl/X1s-exp_i1;
    DI1=const_HI_A2s_1s*(DI1+DnLj);
#endif
    
    if(show_messages>=2) cout << " 2s-1s-2g " << z << " DI1= " << DI1 << endl;
    if(output) Ifile << z << " " << DI1 << endl;
    
    //==========================================================================
    // save
    //==========================================================================
    DI1_2s_vec.push_back(DI1);
    
    return;
}

//==================================================================================================
//==================================================================================================
