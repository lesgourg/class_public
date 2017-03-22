//====================================================================================================================
// accuracy for integrator
//====================================================================================================================
double epsrel_HeI=1.0e-5;

//====================================================================================================================
// to compute the integrals over the solution for the photon field
//====================================================================================================================
bool DP_integrals_splines_are_set=0;
int Spline_for_nx_integral_HeI_memindex;

//====================================================================================================================
// integrand for escape integral
//====================================================================================================================
double dnbar_dx_HeI(double x)
{ return calc_spline_JC(x, Spline_for_nx_integral_HeI_memindex); }

//====================================================================================================================
double calc_DP(double xmin, double xmax, double epsabs=1.0e-50)
{
    //================================================================================================================
    // Integrals
    //================================================================================================================
    double epsrel=epsrel_HeI;
    double r=Integrate_using_Patterson_adaptive(xmin, xmax, epsrel, epsabs, dnbar_dx_HeI);
    
    return r;   
}

//====================================================================================================================
double compute_integral_over_several_resonances(int nres, double nu21, vector<double> &xarr, vector<double> &Fwork)
{
    int nnu_pts=xarr.size();

    double xmin=xarr[0];
    double xmax=xarr[nnu_pts-1];

    double a=xmin, b;
    double dx_x=1.0e-4;
    double r=0.0, r1;
    double epsabs=1.0e-50;
    
    //================================================================================
    // setup splines for interpolation function
    //================================================================================
    if(!DP_integrals_splines_are_set) 
    {
        Spline_for_nx_integral_HeI_memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], &Fwork[0], "DR1sj-integral");
        DP_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral_HeI_memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //================================================================================
    // lower boundary may contain large phase space density (better to separate it)
    //================================================================================
    b=a*(1.0+0.001);
    r1=calc_DP(a, b, epsabs);
    r+=r1; epsabs=fabs(r)*epsrel_HeI;
    a=b;
    
    //================================================================================
    // do resonance by resonance, since the photon field is narrow there at late times
    //================================================================================
    for(int k=0; k<nres-1; k++)
    {
        double xres=phi_i_nui1s_HeI(k)/nu21;
        //============================================================================
        // below resonance 
        //============================================================================
        b=xres*(1.0-dx_x);
        r1=calc_DP(a, b, epsabs);
        r+=r1; epsabs=fabs(r)*epsrel_HeI;
        a=b;
        
        //============================================================================
        // across resonance
        //============================================================================
        b=xres*(1.0+dx_x);
        r1=calc_DP(a, b, epsabs);
        r+=r1; epsabs=fabs(r)*epsrel_HeI;
        a=b;
    }
    
    b=xmax*(1.0-dx_x);
    r1=calc_DP(a, b, epsabs);
    r+=r1; epsabs=fabs(r)*epsrel_HeI;
    
    r1=calc_DP(b, xmax, epsabs);
    r+=r1;

    return r;
}

//====================================================================================================================
// corrections to the escape probabilities of the lines
//====================================================================================================================
double compute_DP_HeI(int line, double z, double nu21, Cosmos &cos, Gas_of_HeI_Atoms &HeIA, 
                      vector<double> &xarr, vector<double> &yarr, vector<double> &Fwork, 
                      vector<double> &exp_x, PDE_solver_functions_HeI &PDE_funcs)
{
    int nnu_pts=xarr.size();
    int ires=PDE_funcs.res_index_of_lines[line];
    int i_HeI=get_HeI_index(ires);
    //
    double nui1=phi_i_nui1s_HeI(line);
    double Tg=cos.TCMB(z), Tm=Tg*calc_HI_rho(z);
    //
    double w=HeIA.Get_gw(i_HeI)/HeIA.Get_gw(0);
    double exp_i1=exp(-const_h_kb*nui1/Tg);
    double Xnl=calc_HeI_Xi(z, ires);
    double X1s=calc_HeI_X1s(z);
    //
    double DnLj=Xnl/X1s/w-exp_i1;
    double Dnem=PDE_funcs.HeI_Dnem_i(z, ires);
    
    //================================================================================
    // tabulate Dn(nu) for splines
    //================================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_i1*(1.0/exp_x[k]-1.0); // exp(-xi1)/npl

        Fwork[k]=(nu21*Dnem-yarr[k]*f_x)*phi_i_HeI(line, k)/phi_i_Dnu_HeI(line); 
    }       
    
    //================================================================================
    // carry out integral
    //================================================================================
    double r=compute_integral_over_several_resonances(PDE_funcs.res_index_of_lines.size(), nu21, xarr, Fwork);
    
    //================================================================================
    // convert to escape probabilities
    //================================================================================
    double tau_S=0.5*w*phi_i_Ai1s_HeI(line)*pow(const_cl/nui1, 3)*X1s*cos.NH(z)/FOURPI/cos.H(z); 
    double Peff=r/Dnem, pd=phi_i_pd_HeI(line), p_sc=1.0-pd;
    double PS=p_ij(tau_S);                               // this is the total Sobolev escape probability 
    double taud=tau_S*pd;
    double Pd_t=p_ij(taud);                              // total escape probability Pd
    //
    Voigtprofile_Dawson &p=get_HeI_profile(line);
    double xmin=xarr[0], xmax=xarr[nnu_pts-1];
    double xDmax=p.nu2x(xmax*nu21, Tm),  chimax=p.chi_Int(xDmax, phi_i_aV_HeI(line));
    double xDmin=p.nu2x(xmin*nu21, Tm),  chimin=p.chi_Int(xDmin, phi_i_aV_HeI(line));   
    //
    double Pd_part=(exp(-chimax*taud)-exp(-chimin*taud))/taud; // this is escape probability inside
                                                               // the computational domain.
    Peff+=Pd_t-Pd_part;                                        // correct for difference to total Pd;
                                                               // At low z this can be a ~0.5 % difference.
    double Pesc=pd*Peff/(1.0-p_sc*Peff);
    
    if(show_messages_HeI>=2) cout << " line = " << line << " z= " << z << " DR = " << r 
                                  << " Pesc = " << Pesc << " PS= " << PS 
                                  << " DPesc= " << Pesc-PS << " pd= " << pd << " Pd_eff= " << Peff
                                  << " tau_S= " << tau_S << " Dnem/DnL= " << Dnem/DnLj << endl;
    
    return Pesc-PS;
}

//====================================================================================================================
// compute convolution of Dn_x with HI continuum opacity.
// This gives the change in the net ionization rate by HeI photons.
//====================================================================================================================
double compute_DR_1sc_absorption_HI(double z, double nu21, Cosmos &cos, Gas_of_HeI_Atoms &HeIA, 
                                    vector<double> &xarr, vector<double> &yarr, vector<double> &Fwork, 
                                    vector<double> sig_Lyc_x, PDE_solver_functions_HeI &PDE_funcs)
{
    int nnu_pts=xarr.size();

    //================================================================================
    // tabulate DN(nu) * sigma^Ly-c_nu for splines
    //================================================================================
    double fac=2.0*FOURPI*pow(nu21/const_cl, 2);
    
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++) 
        Fwork[k]=fac*yarr[k]*pow(xarr[k], 2)*sig_Lyc_x[k]; 
    
    //================================================================================
    // carry out integral
    //================================================================================
    double r=compute_integral_over_several_resonances(PDE_funcs.res_index_of_lines.size(), nu21, xarr, Fwork);
    
    if(show_messages_HeI>=2) cout << " compute_DR_1sc_absorption_HI :: z= " << z 
                                  << " DR = " << r << endl;

    return r;
}

//====================================================================================================================
//
// DI1 output function; correction to the HeI singlet 2s-1s channel
//
//====================================================================================================================
double calc_DI1_2s_Patt_HeI(double xmin, double xmax)
{
    //=========================================================================
    // 1s-->2s Integral; only high nu part (nu/nu21>=0.5) matters
    //=========================================================================
    double r=0.0, epsrel=epsrel_HeI, epsabs=1.0e-60;
    
    r=Integrate_using_Patterson_adaptive(xmin, xmax, epsrel, epsabs, dnbar_dx_HeI); 
    
    return r/const_HeI_A2s_1s;   
}

//====================================================================================================================
void compute_DI1_2s_and_dump_it_HeI(double z, double nu21, vector<double> &xarr, vector<double> &yarr, 
                                    Gas_of_HeI_Atoms &HeIA, Cosmos &cos, 
                                    ofstream &Ifile, vector<double> &DI1_2s_vec, 
                                    vector<double> &Fwork, vector<double> &exp_x,
                                    double (*phi_2g)(int k), PDE_solver_functions_HeI &PDE_funcs,
                                    bool output=0)
{
    //==========================================================================
    int nnu_pts=PDE_funcs.emission_2s1s;
    double nu2s1s=HeIA.Sing.Level(2, 0).Get_Dnu_1s2();
    double Tg=cos.TCMB(z);
    double exp_2s1s=exp(-const_h_kb*nu2s1s/Tg);
    double Dnem=PDE_funcs.HeI_Dnem_i(z, 0);
    
    //==========================================================================
    // tabulate Dn(nu) for splines
    //==========================================================================
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=0; k<nnu_pts; k++)
    { 
        double f_x=exp_2s1s*(1.0/exp_x[k]-1.0); 
        
        if(compute_full_2s_1s_corr_HeI) Fwork[k]=(yarr[k]*f_x-nu21*Dnem)*phi_2g(k); 
        else Fwork[k]=-nu21*Dnem*phi_2g(k); 
    }       
    
    //==========================================================================
    // setup splines
    //==========================================================================
    if(!DP_integrals_splines_are_set) 
    {
        Spline_for_nx_integral_HeI_memindex=calc_spline_coeffies_JC(nnu_pts, &xarr[0], &Fwork[0], "DI1-integral");
        DP_integrals_splines_are_set=1;
    }
    else update_spline_coeffies_JC(Spline_for_nx_integral_HeI_memindex, nnu_pts, &xarr[0], &Fwork[0]);
    
    //==========================================================================
    // compute integrals
    //==========================================================================
    double DI1=calc_DI1_2s_Patt_HeI(0.5*nu2s1s/nu21, nu2s1s/nu21);
    DI1=DI1/Dnem+1.0;
    
    if(show_messages_HeI>=2) cout << " HeI 2s-1s-2g " << z << " DI1= " << DI1 << endl;
    if(output) Ifile << z << " " << DI1 << endl;
    
    //==========================================================================
    // save
    //==========================================================================
    DI1_2s_vec.push_back(DI1);
    
    return;
}

//====================================================================================================================
//====================================================================================================================
