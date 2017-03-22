//==============================================================================
// Routines to evaluate the lhs of the coupled system of ODEs for the
// neutral helium recombination. The system is written in the form
// dXi/dt == g(t, Xi), with Xi=N_i/NHtot
//
// Author: Jens Chluba 
// last modification: May 2011 
//==============================================================================
//
#include "ODEdef_Xi.HeI.h"
#include "Sobolev.h"
#include "HeI_Atom.h"
#include "Photoionization_cross_section.h"
#include "Cosmos.h"
#include "physical_consts.h"
#include "routines.h"
#include "Pesc.HI.h"

//==============================================================================
// equivalent of the Ly-a line
//==============================================================================
double Recombination_calc_p_d(double Tg, double A21L, double Ric, Gas_of_HeI_Atoms &HeIA)
{
    double r=Ric;
    
    int nl=HeIA.Get_total_number_of_Levels();   // Sing + Trip
    int n, l, s;
    int np, lp, sp, jp;
    double gw, gwp;
    double A21, efac;
    
    // start with singlet 3s > 2p   
    for(int i=3; i<nl; i++)
    {
        n =HeIA.Get_n(i);
        l =HeIA.Get_l(i);
        s =HeIA.Get_S(i);
        gw=HeIA.Get_gw(i);
        
        //=========================================================
        // down transititions...
        //=========================================================
        for(int m=0; m<HeIA.Get_n_down(i); m++)
        {
            // here level i is the 'upper' level
            np =HeIA.Get_Trans_Data(i ,m).np;
            lp =HeIA.Get_Trans_Data(i ,m).lp;
            sp =HeIA.Get_Trans_Data(i ,m).sp;
            jp =HeIA.Get_Trans_Data(i ,m).jp;
            
            if(np==2 && lp==1 && sp==0)
            {
                gwp=HeIA.Get_Trans_Data(i ,m).gwp;
                //
                A21=HeIA.Get_Trans_Data(i ,m).A21;
                efac=exp_nu(HeIA.Get_Trans_Data(i ,m).Dnu, Tg);
                r+=A21*gw/gwp/(efac-1.0);
                break;
            }
        } 
    }
    
    //=========================================================
    // also add all other possible ways out (to lower levels)
    //=========================================================
    int index=HeIA.Get_Level_index(2, 1, 0, 1);
    for(int m=0; m<HeIA.Get_n_down(index); m++)
    {
        // exclude 2p->1s transition
        if(HeIA.Get_Trans_Data(index, m).np!=1)
        { 
            A21=HeIA.Get_Trans_Data(index, m).A21;
            efac=exp_nu(HeIA.Get_Trans_Data(index, m).Dnu, Tg);  
            r+=A21*(1.0+1.0/(efac-1.0));
        }
    }
    
    return r/A21L/(1.0+r/A21L);
}

double Recombination_calc_p_sc(double Tg, double A21L, double Ric, Gas_of_HeI_Atoms &HeIA)
{ return 1.0-Recombination_calc_p_d(Tg, A21L, Ric, HeIA); }


//==============================================================================
// death-probability for Singlet Ly-n lines
//==============================================================================
double Recombination_calc_p_d_Singlet(int ni, int li, double Tg, 
                                      double A21L, double Ric, Gas_of_HeI_Atoms &HeIA)
{
    double r=Ric;
    double A21, efac;
    
    //=========================================================
    // also add all other possible ways out (to lower levels)
    //=========================================================
    int index=HeIA.Get_Level_index(ni, li, 0, 1);
    for(int m=0; m<HeIA.Get_n_down(index); m++)
    {
        // exclude np->1s transition
        if(HeIA.Get_Trans_Data(index, m).np!=1)
        { 
            A21=HeIA.Get_Trans_Data(index, m).A21;
            efac=exp_nu(HeIA.Get_Trans_Data(index, m).Dnu, Tg);  
            r+=A21*(1.0+1.0/(efac-1.0));
        }
    }

    int nl=HeIA.Get_total_number_of_Levels();   // Sing + Trip
    int n, l, s;
    int np, lp, sp, jp;
    double gw, gwp;
    
    // start with singlet 2s level
    for(int i=1; i<nl; i++)
    {
        n =HeIA.Get_n(i);
        l =HeIA.Get_l(i);
        s =HeIA.Get_S(i);
        gw=HeIA.Get_gw(i);
        
        //=========================================================
        // down transititions from upper level (n, l, s, j)
        //=========================================================
        for(int m=0; m<HeIA.Get_n_down(i); m++)
        {
            // here level i is the 'upper' level
            np =HeIA.Get_Trans_Data(i ,m).np;
            lp =HeIA.Get_Trans_Data(i ,m).lp;
            sp =HeIA.Get_Trans_Data(i ,m).sp;
            jp =HeIA.Get_Trans_Data(i ,m).jp;
            
            if(np==ni && lp==li && sp==0 && jp==1)
            {
                gwp=HeIA.Get_Trans_Data(i ,m).gwp;
                //
                A21=HeIA.Get_Trans_Data(i ,m).A21;
                efac=exp_nu(HeIA.Get_Trans_Data(i ,m).Dnu, Tg);
                r+=A21*gw/gwp/(efac-1.0);
                break;
            }
        } 
    }
    
    return r/A21L/(1.0+r/A21L);
}

double Recombination_calc_p_sc_Singlet(int ni, int li, double Tg, 
                                       double A21L, double Ric, Gas_of_HeI_Atoms &HeIA)
{ return 1.0-Recombination_calc_p_d_Singlet(ni, li, Tg, A21L, Ric, HeIA); }

//==============================================================================
// death-probability for any level
//==============================================================================
double Recombination_calc_p_d_level(int ni, int li, int si, int ji, double Tg, 
                                    double A21L, double Ric, Gas_of_HeI_Atoms &HeIA)
{
    double r=Ric;
    double A21, efac;

    //=========================================================
    // also add all other possible ways out (to lower levels)
    //=========================================================
    int index=HeIA.Get_Level_index(ni, li, si, ji);
    for(int m=0; m<HeIA.Get_n_down(index); m++)
    {
        // exclude np->1s transition
        if(HeIA.Get_Trans_Data(index, m).np!=1)
        { 
            A21=HeIA.Get_Trans_Data(index, m).A21;
            efac=exp_nu(HeIA.Get_Trans_Data(index, m).Dnu, Tg);  
            r+=A21*(1.0+1.0/(efac-1.0));
        }
    }
    
    int nl=HeIA.Get_total_number_of_Levels();   // Sing + Trip
    int n, l, s;
    int np, lp, sp, jp;
    double gw, gwp;
    
    // start with singlet 2s level
    for(int i=1; i<nl; i++)
    {
        n =HeIA.Get_n(i);
        l =HeIA.Get_l(i);
        s =HeIA.Get_S(i);
        gw=HeIA.Get_gw(i);
        
        //=========================================================
        // down transititions from upper level (n, l, s, j)
        //=========================================================
        for(int m=0; m<HeIA.Get_n_down(i); m++)
        {
            // here level i is the 'upper' level
            np =HeIA.Get_Trans_Data(i ,m).np;
            lp =HeIA.Get_Trans_Data(i ,m).lp;
            sp =HeIA.Get_Trans_Data(i ,m).sp;
            jp =HeIA.Get_Trans_Data(i ,m).jp;
            
            if(np==ni && lp==li && sp==si && jp==ji)
            {
                gwp=HeIA.Get_Trans_Data(i ,m).gwp;
                //
                A21=HeIA.Get_Trans_Data(i ,m).A21;
                efac=exp_nu(HeIA.Get_Trans_Data(i ,m).Dnu, Tg);
                r+=A21*gw/gwp/(efac-1.0);
                break;
            }
        } 
    }
    
    return r/A21L/(1.0+r/A21L);
}

double Recombination_calc_p_sc_level(int ni, int li, int si, int ji, double Tg, 
                                     double A21L, double Ric, Gas_of_HeI_Atoms &HeIA)
{ return 1.0-Recombination_calc_p_d_level(ni, li, si, ji, Tg, A21L, Ric, HeIA); }


//==================================================================================
// to read the diffusion correction function for previous calculation from a file
//==================================================================================
static int fcorr_is_loaded=0;
struct fcorr_spline_data
{
    int nz;
    int memindex;
    double zmin, zmax;
};

fcorr_spline_data fcorr_spline;

//==================================================================================
void Load_fcorr(string fname, fcorr_spline_data &D)
{
    cout << " Load_fcorr:: Loading correction factor for DPesc due to HI :\n " << fname << endl;
    ifstream ifile(fname.c_str());
    if(ifile.is_open()!=1){ cout << " Error opening file. Exiting." << endl; exit(0); }
    
    vector<vector<double> > X_Data;
    
    vector<double> vdum(2);
    while(ifile)
    {
        ifile >> vdum[0];
        ifile >> vdum[1];
        
        X_Data.push_back(vdum);
    }
    
    cout << " Load_fcorr:: Number of redshift points: " << X_Data.size()-1 << endl << endl;
    ifile.close();
    
    //allocate memory for splines;
    D.nz=X_Data.size()-1;
    
    double *za=new double[D.nz];
    double *ya=new double[D.nz];
    
    // check ordering
    if(X_Data[0][0]>X_Data[1][0])
    {
        // now compute all the spline arrays
        for(int i=0; i<D.nz; i++) za[D.nz-1-i]=X_Data[i][0];
        // fcorr
        for(int i=0; i<D.nz; i++) ya[D.nz-1-i]=X_Data[i][1];
    }
    else
    {
        // now compute all the spline arrays
        for(int i=0; i<D.nz; i++) za[i]=X_Data[i][0];
        // fcorr
        for(int i=0; i<D.nz; i++) ya[i]=X_Data[i][1];
    }
    D.zmin=min(za[0], za[D.nz-1]);
    D.zmax=max(za[0], za[D.nz-1]);
    
    D.memindex=calc_spline_coeffies_JC(D.nz, za, ya, "Load_fcorr: fcorr");
    
    delete [] za;
    delete [] ya;
    
    X_Data.clear();
    
    return;      
}

//==================================================================================
double fcorr_Loaded(double Tg)
{
    double z=Tg/2.725-1.0;
    
    if(fcorr_is_loaded==0)
    { 
        Load_fcorr(COSMORECDIR+"./Development/Recombination/Data.fcorr/f.corr.dat", fcorr_spline); 
        fcorr_is_loaded=1; 
    }
    
    if(z>fcorr_spline.zmax) z=fcorr_spline.zmax;
    if(z<fcorr_spline.zmin) z=fcorr_spline.zmin;
    
    return calc_spline_JC(z, fcorr_spline.memindex, "fcorr_Loaded");
}

//==================================================================================
// approximations for DPesc for coherent scattering (no redistribution)
//==================================================================================
double DPesc_coh(double eta_S, double eta_c, 
                 double Tg, double Te, 
                 double A21L, double Ric, Gas_of_HeI_Atoms &HeIA, 
                 Photoionization_cross_section_SH &H1s)
{
    double pd=Recombination_calc_p_d(Tg, A21L, Ric, HeIA);
    double PS=p_ij(eta_S);
    double P_d=p_ij(pd*eta_S);
    double Dpij=DPesc_appr_I_sym(Te, pd*eta_S, eta_c, HeIA.nP_S_profile(2), H1s, P_d*1.0e-5);
    double P=P_d+Dpij;
    double Pesc=pd*P/(1.0-(1.0-pd)*P);
    
    return Pesc-PS;
}

double DPesc_coh(double eta_S, double eta_c, 
                 double Tg, double Te, 
                 double pd, Voigtprofile_Dawson &Prof,
                 Photoionization_cross_section_SH &H1s)
{
    double PS=p_ij(eta_S);
    double P_d=p_ij(pd*eta_S);
    double Dpij=DPesc_appr_I_sym(Te, pd*eta_S, eta_c, Prof, H1s, P_d*1.0e-5);
//  double Dpij=DPesc_appr_I_sym_lg(Te, pd*eta_S, eta_c, Prof, H1s, P_d*1.0e-5);
    double P=P_d+Dpij;
    double Pesc=pd*P/(1.0-(1.0-pd)*P);
    
    return Pesc-PS;
}

//==================================================================================
double DPd_coh(double eta_d, double eta_c, 
               double Tg, double Te, 
               Voigtprofile_Dawson &Prof,
               Photoionization_cross_section_SH &H1s)
{ return DPesc_appr_I_sym(Te, eta_d, eta_c, Prof, H1s, p_ij(eta_d)*1.0e-5); }

//==================================================================================
// using the results obtained with the diffusion code
//==================================================================================
double DPesc_coh_fcorr(double eta_S, double eta_c, 
                       double Tg, double Te, 
                       double A21L, double Ric, Gas_of_HeI_Atoms &HeIA, 
                       Photoionization_cross_section_SH &H1s)
{
    double pd=Recombination_calc_p_d(Tg, A21L, Ric, HeIA);
    double PS=p_ij(eta_S);
    double P_d=p_ij(pd*eta_S);
    double Dpij=DPesc_appr_I_sym(Te, pd*eta_S, eta_c, HeIA.nP_S_profile(2), H1s, P_d*1.0e-5);
    double P=P_d+Dpij;
    double Pesc=pd*P/(1.0-(1.0-pd)*P);
    double fc=fcorr_Loaded(Tg);
    Pesc*=fc;
    
    //cout << z << " " << pd << " " << P << " " << Pesc << " " << Pesc-PS << " " << Pesc/PS-1.0 << endl;
    //cout << z << " fcorr= " << fc << " " << eta_c << " " << eta_S << " " << Dpij << endl;
    
    return Pesc-PS;
}

double DPesc_coh_fcorr(double eta_S, double eta_c, 
                       double Tg, double Te,
                       double pd, Voigtprofile_Dawson &Prof,
                       Photoionization_cross_section_SH &H1s)
{
    double PS=p_ij(eta_S);
    double P_d=p_ij(pd*eta_S);
    double Dpij=DPesc_appr_I_sym(Te, pd*eta_S, eta_c, Prof, H1s, P_d*1.0e-5);
    double P=P_d+Dpij;
    double Pesc=pd*P/(1.0-(1.0-pd)*P);
    double fc=fcorr_Loaded(Tg);
    Pesc*=fc;
    
    //cout << z << " " << pd << " " << P << " " << Pesc << " " << Pesc-PS << " " << Pesc/PS-1.0 << endl;
    //cout << z << " fcorr= " << fc << " " << eta_c << " " << eta_S << " " << Dpij << endl;
    
    return Pesc-PS;
}

//==================================================================================
// to include the effect of HI absorption on the helium 2P - 1S singlet line
//==================================================================================
void evaluate_HI_abs_HeI(const double &z, const double &XH1s, 
                         const double &NH, const double &H_z, 
                         const double &Tg, const double &Tm, 
                         Gas_of_HeI_Atoms &HeI_Atoms, Cosmos &cosmos, 
                         Photoionization_cross_section_SH &H1s, 
                         double (*DPesc_func_ptr)(double Tm, double tauS, double eta_c),
                         double &dXe_HeI_HI_abs_dt)
{
    //------------------------------------------------------------------------------------------
    // need the object HeI_Atoms to contain all the population fractions
    // Xi=NHeI_i/NH at redshift z in addition the photoionization rates
    // have to be up to date.
    //------------------------------------------------------------------------------------------
    // dXi_HeI_dt: will contain all derivatives for the Singlet-states
    // followed by those for the Triplet-states
    //------------------------------------------------------------------------------------------
    double Dpij=0.0, tauS;
    double Xi, Xj;
    //
    double gw=3.0, gwp=1.0;
    double A21, efac, dum, eta_c;
    double sig_c=H1s.sig_phot_ion_nuc();
    
    // HeI-2P-Singlet
    Xi=HeI_Atoms.Xi(2);
    // 1s
    Xj =HeI_Atoms.Xi(0);
    A21=HeI_Atoms.Get_Trans_Data(2, 0).A21;
    efac=exp_nu(HeI_Atoms.Get_Trans_Data(2, 0).Dnu, Tg);
    
    dum=Xi*efac-Xj*gw/gwp;
    //
    tauS=tau_S_gw(A21, HeI_Atoms.Get_Trans_Data(2, 0).lambda21, gw, Xi*NH, gwp, Xj*NH, H_z);
    eta_c=const_cl*NH*XH1s*sig_c/H_z;
    //  
    Dpij=DPesc_func_ptr(Tm, tauS, eta_c);
    
    dXe_HeI_HI_abs_dt=Dpij*A21/(efac-1.0)*dum;
    
//    ofstream ofile;
//    ofile.open("./temp/DP.Sing.dat", ios_base::app);
//    ofile << z << " " << p_ij(tauS)+Dpij << " " << p_ij(tauS) << " " << Dpij << endl;
//    ofile.close();
    
    return;
}


//==================================================================================
// to include the effect of HI absorption on the helium  2P - 1S triplet line
//==================================================================================
void evaluate_HI_abs_HeI_Intercombination(const double &z, const double &XH1s, 
                                          const double &NH, const double &H_z, 
                                          const double &Tg, const double &Tm, 
                                          Gas_of_HeI_Atoms &HeI_Atoms, Cosmos &cosmos, 
                                          Photoionization_cross_section_SH &H1s, 
                                          double (*DPesc_func_ptr)(double Tm, double tauS, double eta_c),
                                          double &dXe_HeI_HI_abs_dt)
{
    //------------------------------------------------------------------------------------------
    // need the object HeI_Atoms to contain all the population fractions
    // Xi=NHeI_i/NH at redshift z in addition the photoionization rates
    // have to be up to date.
    //------------------------------------------------------------------------------------------
    // dXi_HeI_dt: will contain all derivatives for the Singlet-states
    // followed by those for the Triplet-states
    //------------------------------------------------------------------------------------------
    double Dpij, tauS;
    double Xi, Xj;
    //
    double gw=3.0, gwp=1.0;
    unsigned int iT=HeI_Atoms.Get_indexT();
    double A21, efac, dum, eta_c;
    double sig_c=H1s.sig_phot_ion_nuc();
    
    // HeI 2^3 P
    Xi=HeI_Atoms.Xi(2, 1, 1, 1);
    // 1S
    Xj =HeI_Atoms.Xi(0);
    A21=HeI_Atoms.Get_Trans_Data(iT+2, 1).A21;
    efac=exp_nu(HeI_Atoms.Get_Trans_Data(iT+2, 1).Dnu, Tg);
    
    dum=Xi*efac-Xj*gw/gwp;
    //
    tauS=tau_S_gw(A21, HeI_Atoms.Get_Trans_Data(iT+2, 1).lambda21, gw, Xi*NH, gwp, Xj*NH, H_z);
    eta_c=const_cl*NH*XH1s*sig_c/H_z;
    //
    Dpij=DPesc_func_ptr(Tm, tauS, eta_c);

    dXe_HeI_HI_abs_dt=Dpij*A21/(efac-1.0)*dum;
    
//    ofstream ofile;
//    ofile.open("./temp/DP.Trip.dat", ios_base::app);
//    ofile << z << " " << p_ij(tauS)+Dpij << " " << p_ij(tauS) << " " << Dpij << endl;
//    ofile.close();

    return;
}

//------------------------------------------------------------------------------------------
// feedback of Helium lines
//------------------------------------------------------------------------------------------
void evaluate_HeI_feedback(double z, double XH1s, Gas_of_HeI_Atoms &HeI_Atoms, 
                           Cosmos &cosmos, vector<feedback_data > &L,
                           double (*Dn_func)(int, double, double, Gas_of_HeI_Atoms&), 
                           double *dXi_HeI_dt)
{
    //--------------------------------------------------------------------------------------
    // need the object HeI_Atoms to contain all the population fractions
    // Xi=NHeI_i/NH at redshift z
    //--------------------------------------------------------------------------------------
    double NH=cosmos.NH(z), H_z=cosmos.H(z);
    double pij, DRij, tauS, Dn;

    //=========================================================
    // transititions
    //=========================================================
    // update the other way around
    for(int m=(int)L.size()-1; m>=1; m--)
    {
        // evaluate Dn due to distortion from higher line
        Dn=Dn_func(m, z, XH1s, HeI_Atoms);
        
        tauS=tau_S_gw(L[m].A21, L[m].lambda21, 
                      L[m].upperLevel.gw, HeI_Atoms.Xi(L[m].upperLevel.i)*NH, 
                      L[m].lowerLevel.gw, HeI_Atoms.Xi(L[m].lowerLevel.i)*NH, H_z);
        pij =p_ij(tauS);
        DRij=pij*HeI_Atoms.Xi(L[m].lowerLevel.i)*L[m].A21*L[m].upperLevel.gw/L[m].lowerLevel.gw*Dn;
        
        dXi_HeI_dt[L[m].lowerLevel.i]+=-DRij; // remove electron from lower level
        dXi_HeI_dt[L[m].upperLevel.i]+= DRij; // add it to nP-state
    }

    return;
}


