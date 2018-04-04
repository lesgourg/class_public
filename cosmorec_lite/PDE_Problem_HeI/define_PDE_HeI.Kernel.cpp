//==================================================================================================================================
// routines to account for corrections to electron scattering coming form the Doppler Kernel
//==================================================================================================================================
double dmax_kernel=sqrt(-log(1.0e-12));

//==================================================================================================================================
//
// pointer to chosen kernel
//
//==================================================================================================================================
double (*P_Kernel)(double nu, double nup, double Th_e)=PK_Kernel;
//double (*P_Kernel)(double nu, double nup, double Th_e)=P0_Kernel;
//double (*P_Kernel)(double nu, double nup, double Th_e)=P0_Kernel_AliHaimoud;
//double (*P_Kernel)(double nu, double nup, double Th_e)=P_Compton;

//==================================================================================================================================
//
// compute kernel factors for given ix in xarr; This routine is used to illustrate the Kernel.
//
//==================================================================================================================================
void compute_kernel_factors(int ix, double nu21, double Tm, 
                            const vector<double> &xarr, 
                            vector<double> &Kernel_p, 
                            vector<double> &Kernel_m, 
                            double eps)
{
    double Th_e=const_kb_mec2*Tm;
    double nu=xarr[ix]*nu21, fac=nu21/Th_e;
    //
    int np=xarr.size(), k;
    double P0=fac*P_Kernel(nu, nu, Th_e);
    
    //=========================================================================================
    // element nu -> nu
    //=========================================================================================
    Kernel_p[ix]=Kernel_m[ix]=P0;
    
    //=========================================================================================
    // element nu -> nup for nup > nu
    //=========================================================================================
    for(k=ix+1; k<np; k++)
    {
        double nup=xarr[k]*nu21;
        double Dxe=Kernel_deltaxe(nu, nup, Tm);
        
        Kernel_m[k]=fac*P_Kernel(nu, nup, Th_e);
        Kernel_p[k]=Kernel_m[k]*exp(Dxe);
        
        if(Kernel_p[k]/P0<eps)
        {
            for(; k<np; k++) Kernel_p[k]=Kernel_m[k]=0.0;
            break;
        }
    }
    
    //=========================================================================================
    // element nu -> nup for nup < nu
    //=========================================================================================
    for(k=ix-1; k>=0; k--)
    {
        double nup=xarr[k]*nu21;
        double Dxe=Kernel_deltaxe(nu, nup, Tm);
        
        Kernel_m[k]=fac*P_Kernel(nu, nup, Th_e);
        Kernel_p[k]=Kernel_m[k]*exp(Dxe);
        
        if(Kernel_p[k]/P0<eps)
        {
            for(; k>=0; k--) Kernel_p[k]=Kernel_m[k]=0.0;
            break;
        }
    }

    return;
}


//==================================================================================================================================
//
// compute kernel integral over dx for given ix in xarr; Here polynomial interpolation is used for the integrals
//
//==================================================================================================================================
int Kernel_polint_npol=6; // number of points used in interpolation

struct Parameters_for_Kernel_integral
{
    int np, istart;
    vector<double> xarr;
    vector<double> yarr;
    double nu21, nuc, xc, Tm, Th_e, Dnxc;
};

//==================================================================================================================================
// integrand
//==================================================================================================================================
double dKernel_func_polint(double x, void *p)
{ 
    Parameters_for_Kernel_integral *V=((Parameters_for_Kernel_integral *) p);   
    
    if(x<= V->xarr[0] || x>= V->xarr[V->np-1] ) return 0.0;
    
    double Dnx, Dnx_c;
    polint_JC(&V->xarr[0], &V->yarr[0], V->np, x, V->istart, Kernel_polint_npol, &Dnx, &Dnx_c);
    
    double nup=V->nu21*x;
    double Dxe=Kernel_deltaxe(V->nuc, nup, V->Tm);
    double P=P_Kernel(V->nuc, nup, V->Th_e);

    return P*( exp(Dxe)*Dnx - V->Dnxc ); 
}

//==================================================================================================================================
// flip region nup<nu0 to region nup>=nu0
//==================================================================================================================================
double dKernel_func_polint_sym(double dx, void *p)
{ 
    Parameters_for_Kernel_integral *V=((Parameters_for_Kernel_integral *) p);    
    return dKernel_func_polint(V->xc+dx, p)+dKernel_func_polint(V->xc-dx, p); 
}

//==================================================================================================================================
double compute_kernel_Integral_polint(int ix, double nu21, double Tm, 
                                      const vector<double> &xi, 
                                      const vector<double> &yi)
{
    double Th_e=const_kb_mec2*Tm;
    double fac=nu21/Th_e, dx=sqrt(Th_e);
    //
    int np=xi.size(), ndx=2.0*dmax_kernel, iK=ix; // ~ 20.0
    if(iK==1) iK++;
    if(iK==np-2) iK--;
    
    Parameters_for_Kernel_integral Vars;
    Vars.istart=0;
    Vars.xc  =xi[ix];
    Vars.nuc =xi[ix]*nu21;
    Vars.Dnxc=yi[ix];
    Vars.Tm=Tm;
    Vars.Th_e=Th_e;
    Vars.xarr=xi;
    Vars.yarr=yi;
    Vars.nu21=nu21;
    Vars.np=np;
    void *params=&Vars;
    
    //if(Vars.xc>=0.999) return 0.0;
    
    double epsrel=1.0e-6, epsabs=1.0e-60;  
    
    double dx1=min(xi[np-1]-xi[ix], xi[ix]-xi[0]);
    double dxa=min(dx1, xi[ix]*dx*ndx);
    //
    double r=Integrate_using_Patterson_adaptive(0.0, dxa, epsrel, epsabs, dKernel_func_polint_sym, params); 
    //
    r+=Integrate_using_Patterson_adaptive(dxa, xi[ix]*dx*ndx, epsrel, epsabs, dKernel_func_polint_sym, params); 
    return fac*r;
/*    
    double r=Integrate_using_Patterson_adaptive(0.0, xi[ix]*3*dx, epsrel, epsabs, dKernel_func_polint_sym, params); 
    for(int k=2; 3*k<=ndx; k++)
    {
        epsabs=fabs(r)*epsrel;
        double Dr=Integrate_using_Patterson_adaptive(xi[ix]*3*(k-1)*dx, xi[ix]*3*k*dx, epsrel, epsabs, dKernel_func_polint_sym, params);
        r+=Dr;
        
        if(fabs(Dr/r)<epsrel/10.0) break; // stop integrating the decaying kernel part
    }
    
    return r*fac;
*/ 
}

//==================================================================================================================================
//
// routines to account for corrections to electron scattering using Kernel
//
//==================================================================================================================================
void electron_kernel_corr_func(double z, const vector<double> &xi, const vector<double> &yi, 
                               vector<double> &corr_v)
{
//    for(int k=0; k<(int)xi.size(); k++) corr_v[k]=0.0;
//    return;
//    if(z>=2990) return;

    //----------------------------------------------------------
    // xi = nu/nu_Lya
    //----------------------------------------------------------
    double nu21=PDE_Setup_D_HeI.HeI->nP_S_profile(2).Get_nu21();
    
    //----------------------------------------------------------
    // time-dependent variables 
    //----------------------------------------------------------
    double Hz=PDE_Setup_D_HeI.cosm->H(z);
    double NH=PDE_Setup_D_HeI.cosm->NH(z);
    double Tg=PDE_Setup_D_HeI.cosm->TCMB(z);
    
    //----------------------------------------------------------
    // solutions from rec-code
    //----------------------------------------------------------
    double Ne=NH*PDE_Setup_D_HeI.PDE_funcs->HI_Xe(z);
    double rho=PDE_Setup_D_HeI.PDE_funcs->HI_rho(z);
    double Te=Tg*rho;    
    double kap_e=-const_sigT*Ne*(const_kB*Te/const_me_gr/const_cl)/Hz/(1.0+z);
    
    corr_v[0]=corr_v[xi.size()-1]=0.0;
 
#ifdef OPENMP_ACTIVATED    
#pragma omp parallel for default(shared)	
#endif	
    for(int k=1; k<(int)xi.size()-1; k++) 
        corr_v[k]=kap_e*compute_kernel_Integral_polint(k, nu21, Te, xi, yi);
    
    return;
}

//==================================================================================================================================
//==================================================================================================================================



