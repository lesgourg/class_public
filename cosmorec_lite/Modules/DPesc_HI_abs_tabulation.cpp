//===========================================================================================================
// Module to compute the correction to the HeI 2P-1s Singlet and Triplet lines caused by HI 
// absorption. The simple 1D integral approximation of Rubino-Martin et all 2008 is used.
//===========================================================================================================

struct tauS_Block
{ vector<double > DP; };

struct T_Blocks
{
    vector<double > eta;   // the Singlet && Triplet data are ordered according to eta
    vector<double > tauS_S; // within each Singlet[i] the ordering is then according to tauS_x[j]
    vector<tauS_Block > Singlet;
    vector<double > tauS_T;
    vector<tauS_Block > Triplet;
};

struct DPesc_data
{
    vector<double > T;
    vector<T_Blocks> Data;
};

DPesc_data DP_collection;

void save_DP_Data(string filename=Rec_database_path+"/Pesc_Data/DP_Coll_Data.dat", int show=0);
void read_DP_Data(string filename=Rec_database_path+"/Pesc_Data/DP_Coll_Data.dat", int show=0);

static int show_error_message=0;

//===========================================================================================================
// death probablity for 2p Singlet state; this is computed with the effective rates
//===========================================================================================================
double pd_2P_Singlet_HeI(double Tg)
{
    double Bitot;
    get_Bitot_HeI(Tg, Bitot, get_res_index_HeI(2));    
    return 1.0/(1.0+HeI_Atoms.nP_S_profile(2).Get_A21()/Bitot);
}

double pd_2P_Triplet_HeI(double Tg)
{
    return 1.0; // this does not make any significant difference 
    
    //-------------------------------------------------------------    
    double Bitot;
    get_Bitot_HeI(Tg, Bitot, get_res_index_HeI(HeI_Atoms.Get_Level_index(2, 1, 1, 1)));
    
    return 1.0/(1.0+HeI_Atoms.nP_T_profile(2).Get_A21()/Bitot);
}

//===========================================================================================================
// explicit calls
//===========================================================================================================
double call_DP_Singlet(double T, double eta, double tauS, double pd)
{
    if(show_error_message==1) cerr << " Computing Singlet-integral explicitly! T= " 
                                   << T << " eta= " << eta << " tauS= " << tauS*pd << "\n";
    
    return DPesc_coh(tauS, eta, T, T, pd, HeI_Atoms.nP_S_profile(2), HeI_Atoms.HILyc);
}

double call_DP_Triplet(double T, double eta, double tauS, double pd)
{
    if(show_error_message==1) cerr << " Computing Triplet-integral explicitly! T= " 
                                   << T << " eta= " << eta << " tauS= " << tauS*pd << "\n";

    return DPesc_coh(tauS, eta, T, T, pd, HeI_Atoms.nP_T_profile(2), HeI_Atoms.HILyc);
}

//===========================================================================================================
// Computing DP [Singlet]
//===========================================================================================================
double DP_interpol_S_tau(unsigned long int iT, double T, unsigned long int ieta, 
                         double eta, double tauS, double pd)
{
    unsigned long int ntau=DP_collection.Data[iT].tauS_S.size();
    unsigned long int itau;
    double lgtau=log(tauS*pd);
    
    //===========================================================================
    // find index for tau (tau is ordered tau_0 < tau_1 < ... < tau_n)
    //===========================================================================
    double lgtmin=min(DP_collection.Data[iT].tauS_S[0], DP_collection.Data[iT].tauS_S[ntau-1]), 
           lgtmax=max(DP_collection.Data[iT].tauS_S[0], DP_collection.Data[iT].tauS_S[ntau-1]);
    
    if(lgtau<lgtmin || lgtau>lgtmax)
    { 
        if(show_error_message==1) cerr << " tau not inside table " 
                                       << exp(DP_collection.Data[iT].tauS_S[0]) << " " 
                                       << exp(DP_collection.Data[iT].tauS_S[ntau-1]) << " " 
                                       << tauS*pd << "\n"; 
        
        return call_DP_Singlet(T, eta, tauS, pd);
    }
    
    double Dlgt=DP_collection.Data[iT].tauS_S[1]-DP_collection.Data[iT].tauS_S[0];
    double DD=lgtau-DP_collection.Data[iT].tauS_S[0];
    itau=floor(DD/Dlgt);
    
    //===========================================================================
    // region around index
    //===========================================================================
    if(itau>1) itau--;
    if(itau>ntau-4) itau=ntau-4;
    
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double P0=DP_collection.Data[iT].Singlet[ieta].DP[itau];
    double P1=DP_collection.Data[iT].Singlet[ieta].DP[itau+1];
    double P2=DP_collection.Data[iT].Singlet[ieta].DP[itau+2];
    double P3=DP_collection.Data[iT].Singlet[ieta].DP[itau+3];
    //
    double a[4], DD3=pow(Dlgt, 3);
    //
    double D1=lgtau-DP_collection.Data[iT].tauS_S[itau];
    double D2=lgtau-DP_collection.Data[iT].tauS_S[itau+1];
    double D3=lgtau-DP_collection.Data[iT].tauS_S[itau+2];
    double D4=lgtau-DP_collection.Data[iT].tauS_S[itau+3];
    //
    a[0]= D2*D3*D4/(-6.0*DD3);
    a[1]= D1*D3*D4/( 2.0*DD3);
    a[2]= D1*D2*D4/(-2.0*DD3);     
    a[3]= D1*D2*D3/( 6.0*DD3);
    
    //===========================================================================
    return a[0]*P0+a[1]*P1+a[2]*P2+a[3]*P3;
}

double DP_interpol_S_eta(unsigned long int iT, double T, double eta, double tauS, double pd)
{
    unsigned long int neta=DP_collection.Data[iT].eta.size();
    unsigned long int ieta;
    double lgeta=log(eta);
    
    //===========================================================================
    // find index for eta (eta is ordered eta_0 < eta_1 < ... < eta_n)
    //===========================================================================
    double lgemin=min(DP_collection.Data[iT].eta[0], DP_collection.Data[iT].eta[neta-1]), 
           lgemax=max(DP_collection.Data[iT].eta[0], DP_collection.Data[iT].eta[neta-1]);
    
    if(lgeta<lgemin || lgeta>lgemax)
    { 
        if(show_error_message==1) cerr << " eta not inside table " 
                                       << exp(DP_collection.Data[iT].eta[0]) << " " 
                                       << exp(DP_collection.Data[iT].eta[neta-1]) << " " 
                                       << eta << "\n"; 
        
        return call_DP_Singlet(T, eta, tauS, pd);
    }
    
    double Dlge=DP_collection.Data[iT].eta[1]-DP_collection.Data[iT].eta[0];
    double DD=lgeta-DP_collection.Data[iT].eta[0];
    ieta=floor(DD/Dlge);
    
    //===========================================================================
    // region around index
    //===========================================================================
    if(ieta>1) ieta--;
    if(ieta>neta-4) ieta=neta-4;
    
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double P0=DP_interpol_S_tau(iT, T, ieta,   eta, tauS, pd);
    double P1=DP_interpol_S_tau(iT, T, ieta+1, eta, tauS, pd);
    double P2=DP_interpol_S_tau(iT, T, ieta+2, eta, tauS, pd);
    double P3=DP_interpol_S_tau(iT, T, ieta+3, eta, tauS, pd);
    //
    double a[4], DD3=pow(Dlge, 3);
    //
    double D1=lgeta-DP_collection.Data[iT].eta[ieta];
    double D2=lgeta-DP_collection.Data[iT].eta[ieta+1];
    double D3=lgeta-DP_collection.Data[iT].eta[ieta+2];
    double D4=lgeta-DP_collection.Data[iT].eta[ieta+3];
    //
    a[0]= D2*D3*D4/(-6.0*DD3);
    a[1]= D1*D3*D4/( 2.0*DD3);
    a[2]= D1*D2*D4/(-2.0*DD3);     
    a[3]= D1*D2*D3/( 6.0*DD3);
    
    //===========================================================================
    return a[0]*P0+a[1]*P1+a[2]*P2+a[3]*P3;
}

double DP_interpol_S(double T, double tauS, double eta)
{
//    return call_DP_Singlet(T, eta, tauS, pd_2P_Singlet_HeI(T));
//    return 0.0;
    
    unsigned long int nT=DP_collection.T.size();
    unsigned long int iT;
    double lgT=log(T);
    double pd=pd_2P_Singlet_HeI(T);
    
    //===========================================================================
    // find index for Tg (Tg is ordered Tg_0 > Tg_1 > ... > Tg_n)
    //===========================================================================
    double lgTmin=min(DP_collection.T[0], DP_collection.T[nT-1]), 
           lgTmax=max(DP_collection.T[0], DP_collection.T[nT-1]);
    
    if(lgT<lgTmin || lgT>lgTmax)
    { 
        if(show_error_message==1) cerr << " T not inside table " 
                                       << exp(DP_collection.T[0]) << " " 
                                       << exp(DP_collection.T[nT-1]) << " " 
                                       << T << "\n";        
        
        return call_DP_Singlet(T, fabs(eta), fabs(tauS), pd);
    }
    
    locate_JC(&DP_collection.T[0], nT, lgT, &iT); 
    
    //===========================================================================
    // region around index
    //===========================================================================
    if(iT>1) iT--;
    if(iT>nT-4) iT=nT-4;

    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double DPi  =DP_interpol_S_eta(iT,   T, fabs(eta), fabs(tauS), pd);
    double DPip1=DP_interpol_S_eta(iT+1, T, fabs(eta), fabs(tauS), pd);
    double DPip2=DP_interpol_S_eta(iT+2, T, fabs(eta), fabs(tauS), pd);
    double DPip3=DP_interpol_S_eta(iT+3, T, fabs(eta), fabs(tauS), pd);
    //
    double T1 = DP_collection.T[iT];
    double T2 = DP_collection.T[iT+1];
    double T3 = DP_collection.T[iT+2];
    double T4 = DP_collection.T[iT+3];
    //
    double a[4];
    //
    a[0]= (lgT - T2)*(lgT - T3)*(lgT - T4)/((T1 - T2)*(T1 - T3)*(T1 - T4));
    a[1]= (lgT - T1)*(lgT - T3)*(lgT - T4)/((T1 - T2)*(T2 - T3)*(T4 - T2));
    a[2]= (lgT - T1)*(lgT - T2)*(lgT - T4)/((T1 - T3)*(T2 - T3)*(T3 - T4));     
    a[3]= (lgT - T1)*(lgT - T2)*(lgT - T3)/((T3 - T4)*(T4 - T1)*(T2 - T4));
        
    double fc=(_HI_abs_appr_flag==1 ? fcorr_Loaded(T) : 1.0);   
    double DP=a[0]*DPi+a[1]*DPip1+a[2]*DPip2+a[3]*DPip3;
    //
    double PS=p_ij(tauS);
    double P=p_ij(pd*tauS)+DP;
    double Pesc=fc*pd*P/(1.0-(1.0-pd)*P);
    //
    return Pesc-PS;
}

//===========================================================================================================
// Computing DP [Triplet]
//===========================================================================================================
double DP_interpol_T_tau(unsigned long int iT, double T, unsigned long int ieta, 
                         double eta, double tauS, double pd)
{
    unsigned long int ntau=DP_collection.Data[iT].tauS_T.size();
    unsigned long int itau;
    double lgtau=log(tauS*pd);
    
    //===========================================================================
    // find index for tau (tau is ordered tau_0 < tau_1 < ... < tau_n)
    //===========================================================================
    double lgtmin=min(DP_collection.Data[iT].tauS_T[0], DP_collection.Data[iT].tauS_T[ntau-1]), 
           lgtmax=max(DP_collection.Data[iT].tauS_T[0], DP_collection.Data[iT].tauS_T[ntau-1]);
    
    if(lgtau<lgtmin || lgtau>lgtmax)
    { 
        if(show_error_message==1) cerr << " tau not inside table " 
                                       << exp(DP_collection.Data[iT].tauS_T[0]) << " " 
                                       << exp(DP_collection.Data[iT].tauS_T[ntau-1]) << " " 
                                       << tauS*pd << "\n"; 
        
        return call_DP_Triplet(T, eta, tauS, pd);
    }
    
    double Dlgt=DP_collection.Data[iT].tauS_T[1]-DP_collection.Data[iT].tauS_T[0];
    double DD=lgtau-DP_collection.Data[iT].tauS_T[0];
    itau=floor(DD/Dlgt);
    
    //===========================================================================
    // region around index
    //===========================================================================
    if(itau>1) itau--;
    if(itau>ntau-4) itau=ntau-4;
    
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double P0=DP_collection.Data[iT].Triplet[ieta].DP[itau];
    double P1=DP_collection.Data[iT].Triplet[ieta].DP[itau+1];
    double P2=DP_collection.Data[iT].Triplet[ieta].DP[itau+2];
    double P3=DP_collection.Data[iT].Triplet[ieta].DP[itau+3];
    //
    double a[4], DD3=pow(Dlgt, 3);
    //
    double D1=lgtau-DP_collection.Data[iT].tauS_T[itau];
    double D2=lgtau-DP_collection.Data[iT].tauS_T[itau+1];
    double D3=lgtau-DP_collection.Data[iT].tauS_T[itau+2];
    double D4=lgtau-DP_collection.Data[iT].tauS_T[itau+3];
    //
    a[0]= D2*D3*D4/(-6.0*DD3);
    a[1]= D1*D3*D4/( 2.0*DD3);
    a[2]= D1*D2*D4/(-2.0*DD3);     
    a[3]= D1*D2*D3/( 6.0*DD3);
    
    //===========================================================================
    return a[0]*P0+a[1]*P1+a[2]*P2+a[3]*P3;
}

double DP_interpol_T_eta(unsigned long int iT, double T, double eta, double tauS, double pd)
{
    unsigned long int neta=DP_collection.Data[iT].eta.size();
    unsigned long int ieta;
    double lgeta=log(eta);
    
    //===========================================================================
    // find index for eta (eta is ordered eta_0 < eta_1 < ... < eta_n)
    //===========================================================================
    double lgemin=min(DP_collection.Data[iT].eta[0], DP_collection.Data[iT].eta[neta-1]), 
           lgemax=max(DP_collection.Data[iT].eta[0], DP_collection.Data[iT].eta[neta-1]);
    
    if(lgeta<lgemin || lgeta>lgemax)
    { 
        if(show_error_message==1) cerr << " eta not inside table " 
                                       << exp(DP_collection.Data[iT].eta[0]) << " " 
                                       << exp(DP_collection.Data[iT].eta[neta-1]) << " " 
                                       << eta << "\n";  
        
        return call_DP_Triplet(T, eta, tauS, pd);
    }
    
    double Dlge=DP_collection.Data[iT].eta[1]-DP_collection.Data[iT].eta[0];
    double DD=lgeta-DP_collection.Data[iT].eta[0];
    ieta=floor(DD/Dlge);
    
    //===========================================================================
    // region around index
    //===========================================================================
    if(ieta>1) ieta--;
    if(ieta>neta-4) ieta=neta-4;
    
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double P0=DP_interpol_T_tau(iT, T, ieta,   eta, tauS, pd);
    double P1=DP_interpol_T_tau(iT, T, ieta+1, eta, tauS, pd);
    double P2=DP_interpol_T_tau(iT, T, ieta+2, eta, tauS, pd);
    double P3=DP_interpol_T_tau(iT, T, ieta+3, eta, tauS, pd);
    //
    double a[4], DD3=pow(Dlge, 3);
    //
    double D1=lgeta-DP_collection.Data[iT].eta[ieta];
    double D2=lgeta-DP_collection.Data[iT].eta[ieta+1];
    double D3=lgeta-DP_collection.Data[iT].eta[ieta+2];
    double D4=lgeta-DP_collection.Data[iT].eta[ieta+3];
    //
    a[0]= D2*D3*D4/(-6.0*DD3);
    a[1]= D1*D3*D4/( 2.0*DD3);
    a[2]= D1*D2*D4/(-2.0*DD3);     
    a[3]= D1*D2*D3/( 6.0*DD3);
    
    //===========================================================================
    return a[0]*P0+a[1]*P1+a[2]*P2+a[3]*P3;
}

double DP_interpol_T(double T, double tauS, double eta)
{
//    return call_DP_Triplet(T, fabs(eta), fabs(tauS), pd_2P_Triplet_HeI(T));
//    return 0.0;
    
    unsigned long int nT=DP_collection.T.size();
    unsigned long int iT;
    double lgT=log(T);
    double pd=pd_2P_Triplet_HeI(T);
   
    //===========================================================================
    // find index for Tg (Tg is ordered Tg_0 > Tg_1 > ... > Tg_n)
    //===========================================================================
    double lgTmin=min(DP_collection.T[0], DP_collection.T[nT-1]), 
           lgTmax=max(DP_collection.T[0], DP_collection.T[nT-1]);
    
    if(lgT<lgTmin || lgT>lgTmax)
    { 
        if(show_error_message==1) cerr << " T not inside table " 
                                       << exp(DP_collection.T[0]) << " " 
                                       << exp(DP_collection.T[nT-1]) << " " 
                                       << T << "\n"; 
        
        return call_DP_Triplet(T, fabs(eta), fabs(tauS), pd);
    }
        
    locate_JC(&DP_collection.T[0], nT, lgT, &iT); 

    //===========================================================================
    // region around index
    //===========================================================================
    if(iT>1) iT--;
    if(iT>nT-4) iT=nT-4;
        
    //===========================================================================
    // output with 4 point interpolation
    //===========================================================================
    double DPi  =DP_interpol_T_eta(iT  , T, fabs(eta), fabs(tauS), pd);
    double DPip1=DP_interpol_T_eta(iT+1, T, fabs(eta), fabs(tauS), pd);
    double DPip2=DP_interpol_T_eta(iT+2, T, fabs(eta), fabs(tauS), pd);
    double DPip3=DP_interpol_T_eta(iT+3, T, fabs(eta), fabs(tauS), pd);
    //
    double T1 = DP_collection.T[iT];
    double T2 = DP_collection.T[iT+1];
    double T3 = DP_collection.T[iT+2];
    double T4 = DP_collection.T[iT+3];
    //
    double a[4];
    //
    a[0]= (lgT - T2)*(lgT - T3)*(lgT - T4)/((T1 - T2)*(T1 - T3)*(T1 - T4));
    a[1]= (lgT - T1)*(lgT - T3)*(lgT - T4)/((T1 - T2)*(T2 - T3)*(T4 - T2));
    a[2]= (lgT - T1)*(lgT - T2)*(lgT - T4)/((T1 - T3)*(T2 - T3)*(T3 - T4));     
    a[3]= (lgT - T1)*(lgT - T2)*(lgT - T3)/((T3 - T4)*(T4 - T1)*(T2 - T4));
    
    double DP=a[0]*DPi+a[1]*DPip1+a[2]*DPip2+a[3]*DPip3;
    //
    double PS=p_ij(tauS);
    double P=p_ij(pd*tauS)+DP;
    double Pesc=pd*P/(1.0-(1.0-pd)*P);
    //
    return Pesc-PS;
}

//===========================================================================================================
// read data from previous computation
//===========================================================================================================
vector<vector<double> > XHIHeI1s;

void read_X1s_inf(double z_s, double z_e)
{
  double dum;
  vector<double> vdum(3);
  string filename=Rec_database_path+"/Pesc_Data/HI.HeI.1s.dat";
  ifstream file(filename.c_str());

  while(file)
  {
    for(int i=0; i<3; i++)
    {
      file >> dum;
      vdum[i]=dum;
    }
    if((vdum[0]<=z_s*(1.001)) && (vdum[0]>=z_e/1.001)) XHIHeI1s.push_back(vdum);
  }

  file.close();
  return;
}

//===========================================================================================================
// tabulate DP
//===========================================================================================================
void create_Pesc_tables(double z_s, double z_e, Gas_of_HeI_Atoms &HeIA, Cosmos &cos)
{
    read_X1s_inf(z_s, z_e);
    
    int mess=1;
    int npz=XHIHeI1s.size();
    //
    double z, tauS0, eta0, NH, Hz, Tg, XHI1s, XHeI1s, tau, eta, DP, pd;
    //
    double gw=3, gwp=1;
    double A21=HeIA.Get_Trans_Data(2, 1, 0, 0, 0).A21, 
           lambda21=HeIA.Get_Trans_Data(2, 1, 0, 0, 0).lambda21;
    double sig_c=HeIA.HILyc.sig_phot_ion_nuc();
    //
    double tauS_fac=50.0;
    double eta_fac=50.0;
    //
    int ntau_S=31;
    int ntau_T=31;
    int neta=31;
    int zstep=4;
    //
    vector<double> vtau(ntau_S);
    vector<double> veta(neta);
    //
    T_Blocks zB;
    tauS_Block tB;
    
    DP_collection.T.clear();
    DP_collection.Data.clear();

    if(mess>=1) cout << "\n Tabulating DP for 2^1 P_1 - 1^1 S_0 transition. Using ntS= " 
                     << ntau_S << " neta= " << neta << " points. z-step = " << zstep << endl;

    for(int i=0; i<npz; i+=zstep)
    {
        if(i>=npz) break;
        
        z=XHIHeI1s[i][0];
        XHI1s=XHIHeI1s[i][1];
        XHeI1s=XHIHeI1s[i][2];
        //
        NH=cos.NH(z);
        Hz=cos.H(z);    
        Tg=cos.TCMB(z);
        //
        tauS0=tau_S_gw(A21, lambda21, gw, 1.0e-20, gwp, XHeI1s*NH, Hz);
        pd=pd_2P_Singlet_HeI(Tg);
        tauS0*=pd;
        eta0=const_cl*NH*XHI1s*sig_c/Hz;
        //  
        init_xarr(log(tauS0/tauS_fac), log(tauS0*tauS_fac), &vtau[0], ntau_S, 0, 0);
        init_xarr(log(eta0/eta_fac  ), log(eta0*eta_fac  ), &veta[0], neta  , 0, 0);
        
        DP_collection.T.push_back(log(Tg));
        //
        zB.eta.clear();
        zB.tauS_S.clear();
        zB.Singlet.clear();
        zB.tauS_T.clear();
        zB.Triplet.clear();
        //
        for(int e=0; e<neta; e++) zB.eta.push_back(veta[e]);
        for(int t=0; t<ntau_S; t++) zB.tauS_S.push_back(vtau[t]);
        //
        for(int e=0; e<neta; e++)
        {
            tB.DP.clear();
            eta=exp(veta[e]);
            for(int t=0; t<ntau_S; t++)
            {
                tau=exp(vtau[t]);
                //
                DP=DPd_coh(tau, eta, Tg, Tg, HeI_Atoms.nP_S_profile(2), HeI_Atoms.HILyc);
                tB.DP.push_back(DP);
                
                if(mess>=2) cout << z << " " << eta << " " << tau << " " << DP << endl;
            }
            zB.Singlet.push_back(tB);
        }   
        DP_collection.Data.push_back(zB);       
        
        if(mess==1) cout << "\033[2K" << " z= " << z << " Tg= " << Tg << endl << "\033[1A";
    }

    if(mess>=1) cout << " Tabulating DP for 2^3 P_1 - 1^1 S_0 transition. Using ntT= " 
                     << ntau_T << " neta= " << neta << " points. z-step = " << zstep << endl;

    gw=3; gwp=1;
    A21=HeIA.Get_Trans_Data(HeIA.Get_indexT()+2, 1, 0, 0, 0).A21;
    lambda21=HeIA.Get_Trans_Data(HeIA.Get_indexT()+2, 1, 0, 0, 0).lambda21;
    
    for(int i=0, kl=0; i<npz; i+=zstep, kl++)
    {
        if(i>=npz) break;
        
        z=XHIHeI1s[i][0];
        XHI1s=XHIHeI1s[i][1];
        XHeI1s=XHIHeI1s[i][2];
        //
        NH=cos.NH(z);
        Hz=cos.H(z);    
        Tg=cos.TCMB(z);
        //
        tauS0=tau_S_gw(A21, lambda21, gw, 1.0e-20, gwp, XHeI1s*NH, Hz);
        //  
        vtau.resize(ntau_T);
        init_xarr(log(tauS0/tauS_fac), log(tauS0*tauS_fac), &vtau[0], ntau_T, 0, 0);
        //
        for(int t=0; t<ntau_T; t++) DP_collection.Data[kl].tauS_T.push_back(vtau[t]);
        for(int e=0; e<neta; e++)
        {
            tB.DP.clear();
            eta=exp(DP_collection.Data[kl].eta[e]);
            for(int t=0; t<ntau_T; t++)
            {
                tau=exp(vtau[t]);
                //
                DP=DPd_coh(tau, eta, Tg, Tg, HeI_Atoms.nP_T_profile(2), HeI_Atoms.HILyc);
                tB.DP.push_back(DP);
                
                if(mess>=2) cout << z << " " << eta << " " << tau << " " << DP << endl;
            }
            DP_collection.Data[kl].Triplet.push_back(tB);       
        }   
        if(mess==1) cout << "\033[2K" << " z= " << z << " Tg= " << Tg << endl << "\033[1A";
    }
    
    save_DP_Data(Rec_database_path+"/Pesc_Data/DP_Coll_Data.31.fac_50.neff_"
                                  +int_to_string(nS_effective_HeI)+".dat", mess);

    return;
}

//===========================================================================================================
// Data access
//===========================================================================================================
void save_DP_Data(string filename, int mess)
{
    ofstream file(filename.c_str());
    file.precision(8);
    
    if(mess>=1) cout << " Saving DP-data " << endl;
    
    int nz=DP_collection.T.size(), neta, ntau;
    for(int i=0; i<nz; i++) file << DP_collection.T[i] << " ";
    file << endl;
    for(int i=0; i<nz; i++)
    {
        file << "TB" << endl;
        neta=DP_collection.Data[i].eta.size();
        for(int k=0; k<neta; k++) file << DP_collection.Data[i].eta[k] << " ";
        file << endl << "tS" << endl;       
        ntau=DP_collection.Data[i].tauS_S.size();
        for(int l=0; l<ntau; l++) file << DP_collection.Data[i].tauS_S[l] << " ";
        file << endl << "tT" << endl;       
        ntau=DP_collection.Data[i].tauS_T.size();
        for(int l=0; l<ntau; l++) file << DP_collection.Data[i].tauS_T[l] << " ";
        
        file << endl << "DP_S"<< endl;
        for(int k=0; k<neta; k++) 
        {
            ntau=DP_collection.Data[i].tauS_S.size();
            for(int l=0; l<ntau; l++) file << DP_collection.Data[i].Singlet[k].DP[l] << " ";
            file << endl;
        }
            
        file << "DP_T"<< endl;
        for(int k=0; k<neta; k++) 
        {
            ntau=DP_collection.Data[i].tauS_T.size();
            for(int l=0; l<ntau; l++) file << DP_collection.Data[i].Triplet[k].DP[l] << " ";
            file << endl;
        }
    }
    
    file.close();
    
    return;
}

//===========================================================================================================
vector<double > read_block(ifstream &file, string untilstr, string mess, int show=0)
{
    string str="";
    double val;
    vector<double > vals;

    while(file && str!=untilstr)
    {
        file >> str; val=atof(str.c_str());
        if(str!=untilstr)
        {   
            if(show>=2) cout << str << mess << val << " " << exp(val) << endl;
            vals.push_back(val);
        }
    }
    if(show>=2) cout << str << endl;
    
    return vals;
}

//===========================================================================================================
vector<double > read_block(ifstream &file, int np, string mess, int show=0)
{
    double val;
    vector<double > vals;
    
    for(int i=0; i<np; i++)
    { 
        file >> val;
        if(show>=2) cout << i << mess << val << " " << exp(val) << endl;
        vals.push_back(val);
    }
    
    return vals;
}

//===========================================================================================================
void read_DP_Data(string filename, int show)
{    
    ifstream file(filename.c_str());
    if(!file)
    { 
        cout << " read_DP_Data :: creating DP-file for HI-absorption."
             << " This takes a few minutes but is only done once. " << endl;
        
        create_Pesc_tables(3450.0, 1500.0, HeI_Atoms, cosmos);
        file.open(filename.c_str());
    }
    
    if(show>=1) print_message(" Reading DP-data. This can take a little moment... ");
    
    DP_collection.T.clear();
    DP_collection.Data.clear();
    
    T_Blocks zB;
    tauS_Block tB;
    string str;
    
    DP_collection.T=read_block(file, "TB", " z ", show);
    //
    int nz=DP_collection.T.size();

    zB.eta.clear(); zB.Singlet.clear(); zB.Triplet.clear();
    for(int i=0; i<nz; i++) DP_collection.Data.push_back(zB);
        
    for(int i=0; i<nz; i++)
    {
        DP_collection.Data[i].eta=read_block(file, "tS", " eta ", show);
        DP_collection.Data[i].tauS_S=read_block(file, "tT", " tauS_S ", show);
        DP_collection.Data[i].tauS_T=read_block(file, "DP_S", " tauS_T ", show);

        int neta=DP_collection.Data[i].eta.size();
        int ntau=DP_collection.Data[i].tauS_S.size();
        for(int e=0; e<neta; e++)
        {
            tB.DP=read_block(file, ntau, " DP_S ", show);
            DP_collection.Data[i].Singlet.push_back(tB);
        }       
        file >> str;
        if(show>=2) cout << str << endl;
        if(str!="DP_T"){ cerr << " Error 2: Please check DP-file and i/o \n"; exit(0); }

        ntau=DP_collection.Data[i].tauS_T.size();
        for(int e=0; e<neta; e++)
        {
            tB.DP=read_block(file, ntau, " DP_T ", show);
            DP_collection.Data[i].Triplet.push_back(tB);
        }       
        file >> str;
        if(show>=2) cout << str << endl;
        if(str!="TB" && file){ cerr << " Error 3: Please check DP-file and i/o \n"; exit(0); }
    }
    
    file.close();
        
    return;
}

