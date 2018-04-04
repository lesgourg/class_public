//==================================================================================================
// to include the distortion from HeI-lines
//==================================================================================================
struct Data_Dn
{
    double *z; // redshift
    double *Dn; // Dn_(zem) phase space density at z=zem
    double *zpd; // pd(z) death probability for the level
    double *pd; // pd(z) death probability for the level
};

//==================================================================================================
// feedback all helium lines
//==================================================================================================
vector<feedback_data > HeI_Lines;
int HeI_photons_nz;
int HeI_photons_iznow;
Data_Dn *HeI_photons;

//==================================================================================================
double interpol_HeI_Dn(double x, double *xa, double *ya, int ny)
{  
    if(ny<=4 || x<min(xa[0], xa[ny-1]) || x>max(xa[0], xa[ny-1])) return 0.0;
    
    //======================================================================
    // here arrays are assumed to start at index 0 and to end at ny-1
    //======================================================================
    int npol=4;
    double v, dv;
    polint_JC(xa, ya, ny, x, npol, &v, &dv);
    
    return v;
}

//==================================================================================================
void create_feedback_Table(int nmax, Gas_of_HeI_Atoms &HeIA)
{
    HeI_Lines.clear();
    
    feedback_data dum;
    int inserted;

    //======================================================================
    // empty level info
    //======================================================================
    dum.upperLevel.n=dum.upperLevel.l=dum.upperLevel.s=dum.upperLevel.j=dum.upperLevel.i=0;
    dum.upperLevel.gw=0.0;
    dum.lowerLevel.n=dum.lowerLevel.l=dum.lowerLevel.s=dum.lowerLevel.j=dum.lowerLevel.i=0;
    dum.lowerLevel.gw=0.0;
    dum.A21=dum.lambda21=0.0;
    dum.nu21=HeIA.HILyc.Get_nu_ionization();
    HeI_Lines.push_back(dum);

    //======================================================================
    // all states
    //======================================================================
    for(int lev=0; lev<(int) HeIA.Get_total_number_of_Levels(); lev++)
    {
        if(HeIA.Get_n(lev)>nmax) continue;
        
        for(int k=0; k<HeIA.Get_n_down(lev); k++)
        {
            //======================================================================
            // write down every line that has transition frequency nu>nuLy-aT, i.e. 
            // exclude triplet line
            //======================================================================
            if(5069096.0e+9<HeIA.Get_Trans_Data(lev, k).Dnu)
            {
                dum.upperLevel.n=HeIA.Get_n(lev);
                dum.upperLevel.l=HeIA.Get_l(lev);
                dum.upperLevel.s=HeIA.Get_S(lev);
                dum.upperLevel.j=HeIA.Get_J(lev);
                dum.upperLevel.gw=HeIA.Get_gw(lev);
                dum.upperLevel.i=lev;
                //
                dum.lowerLevel.n=HeIA.Get_Trans_Data(lev, k).np;
                dum.lowerLevel.l=HeIA.Get_Trans_Data(lev, k).lp;
                dum.lowerLevel.s=HeIA.Get_Trans_Data(lev, k).sp;
                dum.lowerLevel.j=HeIA.Get_Trans_Data(lev, k).jp;
                dum.lowerLevel.gw=HeIA.Get_Trans_Data(lev, k).gwp;
                dum.lowerLevel.i=HeIA.Get_Level_index(dum.lowerLevel.n, dum.lowerLevel.l, 
                                                      dum.lowerLevel.s, dum.lowerLevel.j);
                //  
                dum.A21=HeIA.Get_Trans_Data(lev, k).A21;
                dum.lambda21=HeIA.Get_Trans_Data(lev, k).lambda21;
                dum.nu21=HeIA.Get_Trans_Data(lev, k).Dnu;
                
                //===================================================================
                // if there is only one element
                //===================================================================
                if(HeI_Lines.size()==1) HeI_Lines.push_back(dum); 
                //===================================================================
                // otherwise insert where appropriate
                //===================================================================
                else
                {
                    inserted=0;
                    if(dum.nu21<HeI_Lines[1].nu21)
                    {   
                        HeI_Lines.insert(HeI_Lines.begin()+1, dum); 
                        inserted=1; 
                    }
                    else
                    {
                        for(int t=2; t<(int)HeI_Lines.size(); t++)
                        {
                            if(dum.nu21>HeI_Lines[t-1].nu21 && dum.nu21<HeI_Lines[t].nu21)
                            { 
                                HeI_Lines.insert(HeI_Lines.begin()+t, dum); 
                                inserted=1; 
                            }
                        }
                    }
                    
                    if(inserted==0) HeI_Lines.push_back(dum);
                }
            }
        }
    }
    
    cout << " create_feedback_Table:: Helium feedback table. " << endl << endl;

    int k=1;
    cout << " " << k << " " 
         << HeI_Lines[k].lowerLevel.i << " (" 
         << HeI_Lines[k].lowerLevel.n << " " << HeI_Lines[k].lowerLevel.l << " " 
         << HeI_Lines[k].lowerLevel.s << " " << HeI_Lines[k].lowerLevel.j << ") " 
         << HeI_Lines[k].lowerLevel.gw << " " 
         << HeI_Lines[k].upperLevel.i << " (" 
         << HeI_Lines[k].upperLevel.n << " " << HeI_Lines[k].upperLevel.l << " " 
         << HeI_Lines[k].upperLevel.s << " " << HeI_Lines[k].upperLevel.j << ") "
         << HeI_Lines[k].upperLevel.gw << " "
         << HeI_Lines[k].nu21 << " " << HeI_Lines[k].A21 << " Dz/z = -- " << endl;

    //======================================================================
    // other lines 
    //======================================================================
    for(k=2; k<(int)HeI_Lines.size(); k++) 
        cout << " " << k << " " 
             << HeI_Lines[k].lowerLevel.i << " (" 
             << HeI_Lines[k].lowerLevel.n << " " << HeI_Lines[k].lowerLevel.l << " " 
             << HeI_Lines[k].lowerLevel.s << " " << HeI_Lines[k].lowerLevel.j << ") " 
             << HeI_Lines[k].lowerLevel.gw << " " 
             << HeI_Lines[k].upperLevel.i << " (" 
             << HeI_Lines[k].upperLevel.n << " " << HeI_Lines[k].upperLevel.l << " " 
             << HeI_Lines[k].upperLevel.s << " " << HeI_Lines[k].upperLevel.j << ") " 
             << HeI_Lines[k].upperLevel.gw << " "
             << HeI_Lines[k].nu21 << " " << HeI_Lines[k].A21 << " Dz/z = " 
             << HeI_Lines[k-1].nu21/HeI_Lines[k].nu21 - 1.0 << endl;

    cout << endl;
    
    if(_HeI_feedback_w_HI_abs==0) 
        cout << " create_feedback_Table:: HI opacity will NOT be included between the lines " << endl;
    else if(_HeI_feedback_w_HI_abs==1) 
        cout << " create_feedback_Table:: HI opacity will be included between the lines " << endl;

    return;
}

//==================================================================================================
// for other HeI photons
//==================================================================================================
void init_HeI_feedback(Gas_of_HeI_Atoms &HeIA, Parameters_form_initfile &params, double *za)
{
    HeI_photons_iznow=0;
    
    if(!(HeI_Lines.size()>1))
    {
        create_feedback_Table(min(nHeFeedback, 20), HeIA);
        
        cout << " init_HeI_feedback :: setting up pd and memory " << endl;
        
        //===================================================================
        // to include HeI photons
        //===================================================================
        HeI_photons=new Data_Dn[HeI_Lines.size()];
        
        //===================================================================
        // allocate memory
        //===================================================================
        for(int line=0; line<(int) HeI_Lines.size(); line++)
        {
            HeI_photons[line].z=new double[params.nz+5];
            HeI_photons[line].Dn=new double[params.nz+5];
            HeI_photons[line].zpd=new double[params.nz+5];
            HeI_photons[line].pd=new double[params.nz+5];
        }
    }
    
    //===================================================================
    // fill pd-arrays and reset HeI feedback table; This is only done
    // when the cosmology changed.
    //===================================================================
    if(iteration_count==Diff_iteration_min)
    {
        vector<double> Bitot(get_number_of_resolved_levels_HeI());
        for(int k=0; k<params.nz; k++)
        { 
            // Helium is totally recombined at low redshift. Just stop tabulation below.
            if(za[k]<=500.0){ HeI_photons_nz=k-1; break; } 
            
            double Tg=cosmos.TCMB(za[k]);
            get_Bitot_HeI(Tg, Bitot);
            
            //===========================================================
            // compute death probabilities from effective rate tables
            //===========================================================
            for(int line=1; line<(int) HeI_Lines.size(); line++)
            {
                HeI_photons[line].z[k]=za[k];
                HeI_photons[line].zpd[k]=za[k];
                HeI_photons[line].Dn[k]=0.0;            
                HeI_photons[line].pd[k]=log(1.0/(1.0+HeI_Lines[line].A21
                                          /Bitot[get_res_index_HeI(HeI_Lines[line].upperLevel.i)]));          
            }
        }
    }
    
    return;
}

//==================================================================================================
// interpolation function
//==================================================================================================
double Dn_minus_HeI(int line, double zi)
{   // the last point was saves at i=HeI_photons_iznow
    return interpol_HeI_Dn(zi, HeI_photons[line].z, HeI_photons[line].Dn, 
                           min(HeI_photons_nz, max(HeI_photons_iznow+2, 4)));  
}

//==================================================================================================
double pd_HeI(int line, double z)
{ 
    double pd;
    if(z>HeI_photons[line].z[0] || z<HeI_photons[line].z[HeI_photons_nz-1]) 
    { cout << " pd_HeI:: out of array bounds... Exiting. " << endl; exit(0); }  
    else pd=exp(interpol_HeI_Dn(z, HeI_photons[line].zpd, HeI_photons[line].pd, HeI_photons_nz)); 
    
    return pd;
} 

//==================================================================================================
// feedback with tauc
//==================================================================================================
// defined a bit below (recursive)
void update_current_feedback_functions(int line, Gas_of_HeI_Atoms &HeIA, double XHI1s, int iznow, 
                                       double znow, string mess="");

//==================================================================================================
double Saved_XHI1s(double z)
{ 
    int nmax=min(HeI_photons_nz, max(HeI_photons_iznow+2, 4));
    return interpol_HeI_Dn(z, HeI_photons[0].z, HeI_photons[0].Dn, nmax); 
}

//==================================================================================================
// tauc-integral
//==================================================================================================
double FB_HeI_sig_HI_1s(double nu){ return HeI_Atoms.HILyc.sig_phot_ion(nu); }

double FB_HeI_dtau_c_HI1s(double Dz, void *params) 
{
    // x_em=nu/eta(z)
    double *p=((double *) params);
    double x_em=p[0];
    double z=p[1]+Dz;
    
    double nu0=HeI_Atoms.HILyc.Get_nu_ionization();
    double nu=x_em*nu0*(1.0+z);
    
    double Hz=cosmos.H(z);
    double NHI1s=cosmos.NH(z)*Saved_XHI1s(z);   
    double sig_c=FB_HeI_sig_HI_1s(nu);
    
    return sig_c*NHI1s*const_cl/( Hz*(1.0+z) );
}

//==================================================================================================
// Work horse
//==================================================================================================
double FB_HeI_tau_c_HI1s_region_GSL(double (*fptr)(double, void *), double zmin, double zmax, 
                                    double x_em, double z_em, 
                                    double epsabs, double epsrel, int flag, string mess="")
{
    if(zmin>=zmax) return 0.0;
    if(flag!=0) cout << " tauc_HI1s-region: " << mess << " " << zmin << " " << zmax << " " 
                     << x_em << " " << 1.0/x_em-1.0 << " " << z_em << endl;
    
    int neval=50000;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(neval);
    gsl_function F;
    F.function = fptr;
    double para[3]={x_em, zmin, 0.0};
    F.params=para;
    
    double r=0.0;
    double epsest;
    double ab=0.0, bb=zmax-zmin;
    
    gsl_integration_qag (&F, ab, bb, epsabs, epsrel, neval, 2, w, &r, &epsest); 
    gsl_integration_workspace_free(w);
    
    if(fabs(r)<=epsabs) r=0.0;
    if(flag!=0) cout << " tau = " << r << " done " << endl;
    return r;
}

//==================================================================================================
double FB_HeI_tau_c_HI1s_GSL(double y_em, double z_em, double z_obs, double (*fptr)(double, void *))
{
    if(z_em<=z_obs) return 0.0;
    //=========================================================================
    // y_em=nu_em/nu0, i.e. initial frequency at z=z_em
    //=========================================================================
    double x_em=y_em/(1.0+z_em);
    double zc=1.0/x_em-1.0; 
    if(z_em<=zc) return 0.0; // photons emitted below threshold --> tauc=0
    
    double zlow=max(zc, z_obs); // integral only goes from z_em until max(zc, z_obs)
    
    //=========================================================================
    // build intervals
    //=========================================================================
    vector<double> intervals;
    intervals.push_back(zlow);
    intervals.push_back(z_em);
    
    //=========================================================================
    // Integrals
    //=========================================================================
    double r=0.0, r1=0.0;
    double ab, bb;
    double epsrel=1.0e-4, epsabs=1.0e-10;
    for(int i=1; i<(int)intervals.size(); i++)
    {
        ab=intervals[i-1]; bb=intervals[i];
        r1=FB_HeI_tau_c_HI1s_region_GSL(fptr, ab, bb, x_em, z_em, epsabs, epsrel, 
                                        0, "Interval: " + int_to_string(i));
        r+=r1; 
    }
    
    return r;
}

//==================================================================================================
// optical depth in Lyman continuum
// flag==1 mean fnu=1; otherwise fnu!=1 is used
//==================================================================================================
double FB_HeI_tau_c_HI1s(double y_em, double z_em, double z_obs)
{ return FB_HeI_tau_c_HI1s_GSL(y_em, z_em, z_obs, FB_HeI_dtau_c_HI1s); }

//==================================================================================================
double DP_interpol_S(double Tg, double tauS, double eta);  // defined in P_HI_abs_tabulation.cpp
double DP_interpol_T(double Tg, double tauS, double eta);  // defined in P_HI_abs_tabulation.cpp

//==================================================================================================
double Dn_HeI_effective(int linenum, double znow, double X1sHI, 
                        double Dnminus, Gas_of_HeI_Atoms &HeIA)
{
    if(Dnminus==0.0) return 0.0;        
    //=======================================================================
    // evaluate PdS etc
    //=======================================================================
    double NH=cosmos.NH(znow), H_z=cosmos.H(znow), Tg=cosmos.TCMB(znow);
    
    double tauS=tau_S_gw(HeI_Lines[linenum].A21, HeI_Lines[linenum].lambda21, 
                         HeI_Lines[linenum].upperLevel.gw, 
                         HeIA.Xi(HeI_Lines[linenum].upperLevel.i)*NH, 
                         HeI_Lines[linenum].lowerLevel.gw, 
                         HeIA.Xi(HeI_Lines[linenum].lowerLevel.i)*NH, H_z);
    
    double pd=pd_HeI(linenum, znow);
    
    double Pd=p_ij(pd*tauS), PS=p_ij(tauS), psc=1.0-pd;
    double PdS=pd*Pd/(1.0-psc*Pd);
    if(pd*tauS<=1.0e-6) PdS=(1.0-tauS*pd/2.0)/(1.0+psc*tauS/2.0);
    
    //=======================================================================
    // DnL(z)
    //=======================================================================
    double nL=HeI_Lines[linenum].lowerLevel.gw*HeIA.Xi(HeI_Lines[linenum].upperLevel.i)
              /( HeI_Lines[linenum].upperLevel.gw*HeIA.Xi(HeI_Lines[linenum].lowerLevel.i) );
    nL=1.0/(1.0/nL-1.0);
    double nPl=1.0/(exp_nu(HeI_Lines[linenum].nu21, Tg)-1.0);
    double DnL=nL-nPl;
    //
    double Dn_noscatt=(PS-PdS)/PS*DnL;
    double Dn_feedback=PdS/PS*Dnminus;
    
    //=======================================================================
    // include DP correction for feedback
    //=======================================================================
    double DP=0.0;
    double eta_c=const_cl*NH*X1sHI*HeIA.HILyc.sig_phot_ion_nuc()/H_z;

    if(linenum==1 && flag_HI_absorption==1 && znow<=zcrit_HI_Int && flag_spin_forbidden==1) 
        DP=DP_interpol_T(Tg, tauS, eta_c);

    //=======================================================================
    // here tau_S was used in the tabulation NOT pd*tauS!   
    //=======================================================================
    if(linenum==2 && flag_HI_absorption==1 && znow<=zcrit_HI) DP=DP_interpol_S(Tg, tauS, eta_c);

    PdS=PS+DP;  
    Dn_noscatt=0.0;
    Dn_feedback=PdS/PS*Dnminus;

    return Dn_noscatt+Dn_feedback;
}

//==================================================================================================
double Dn_minus_HeI_feedback(int linenum, double z, double X1sHI, Gas_of_HeI_Atoms &HeIA)
{
    if(HeI_photons_iznow<=4) return 0.0;
    if(linenum+1>=(int)HeI_Lines.size() || linenum==0) return 0.0;
    
    double nu21u=HeI_Lines[linenum+1].nu21;
    double nu21l=HeI_Lines[linenum].nu21;
    double zi=nu21u/nu21l*(1.0+z)-1.0;
    double tauc=0.0, Dnminus=0.0;
    
    if(zi<HeI_photons[0].z[1]) // before that no info for feedback is saved...
    {
        //====================================================================
        // write current redshift point to the array
        //====================================================================
        update_current_feedback_functions(linenum+1, HeIA, X1sHI, HeI_photons_iznow+1, z);
        tauc=const_cl*X1sHI*cosmos.NH(z)*HeIA.HILyc.sig_phot_ion(nu21l)/cosmos.H(z)
                     *(1.0-pow(nu21l/nu21u, 3))/3.0;
    
        Dnminus=Dn_minus_HeI(linenum+1, zi)*exp(-tauc);
    }
    
    return Dnminus;
}

double Dn_HeI_feedback(int linenum, double z, double X1sHI, Gas_of_HeI_Atoms &HeIA)
{ return Dn_HeI_effective(linenum, z, X1sHI, Dn_minus_HeI_feedback(linenum, z, X1sHI, HeIA), HeIA); }

//==================================================================================================
// feedback without tauc
//==================================================================================================
double Dn_minus_HeI_feedback_no_HI_abs(int linenum, double z, double X1sHI, Gas_of_HeI_Atoms &HeIA)
{
    if(HeI_photons_iznow<=4) return 0.0;
    if(linenum+1>=(int)HeI_Lines.size() || linenum==0) return 0.0;
    
    double nu21u=HeI_Lines[linenum+1].nu21;
    double nu21l=HeI_Lines[linenum].nu21;
    double zi=nu21u/nu21l*(1.0+z)-1.0;
    double Dnminus=0.0;

    if(zi<HeI_photons[0].z[1]) // before that no info for feedback is saved...
    {
        //====================================================================
        // write current redshift point to the array
        //====================================================================
        update_current_feedback_functions(linenum+1, HeIA, X1sHI, HeI_photons_iznow+1, z);
        
        Dnminus=Dn_minus_HeI(linenum+1, zi);
    }
    
    return Dnminus;
}

double Dn_HeI_feedback_no_HI_abs(int linenum, double z, double X1sHI, Gas_of_HeI_Atoms &HeIA)
{ return Dn_HeI_effective(linenum, z, X1sHI, Dn_minus_HeI_feedback_no_HI_abs(linenum, z, X1sHI, HeIA), HeIA); }

//==================================================================================================
// saving the feedback information
//==================================================================================================
void update_current_feedback_functions(int line, Gas_of_HeI_Atoms &HeIA, double XHI1s, int iznow, 
                                       double znow, string mess)
{
    //====================================================================
    // return to solver if you are just updating the last point or before that
    //====================================================================
    if(znow>=HeI_photons[line].z[iznow-1]) return; 

    //====================================================================
    // Also save redshift
    //====================================================================
    HeI_photons[0].z[iznow]=znow;

    //====================================================================
    // Also save XHI1s
    //====================================================================
    HeI_photons[0].Dn[iznow]=XHI1s;

    //====================================================================
    // in this case only hydrogen was updated
    //====================================================================
    if(line==0) return;
    
    double nu21=HeI_Lines[line].nu21;
    double nL=HeI_Lines[line].lowerLevel.gw*HeIA.Xi(HeI_Lines[line].upperLevel.i)
              /( HeI_Lines[line].upperLevel.gw*HeIA.Xi(HeI_Lines[line].lowerLevel.i) );
    nL=1.0/(1.0/nL-1.0);
    double nPl=1.0/(exp_nu(nu21, cosmos.TCMB(znow))-1.0);
    double DnL=nL-nPl;
        
    //================================================================================
    // also add the distortions from the previous line that may pass through the line
    //================================================================================
    double Dn=0.0, NH=cosmos.NH(znow), H_z=cosmos.H(znow);

    //====================================================================
    // include the absorption in the line
    //====================================================================
    double tauS=tau_S_gw(HeI_Lines[line].A21, HeI_Lines[line].lambda21, 
                         HeI_Lines[line].upperLevel.gw, 
                         HeIA.Xi(HeI_Lines[line].upperLevel.i)*NH, 
                         HeI_Lines[line].lowerLevel.gw, 
                         HeIA.Xi(HeI_Lines[line].lowerLevel.i)*NH, H_z);
    
    double pd=pd_HeI(line, znow);

    if(_HeI_feedback_w_HI_abs==0) Dn=Dn_minus_HeI_feedback_no_HI_abs(line, znow, XHI1s, HeIA);
    else if(_HeI_feedback_w_HI_abs==1) Dn=Dn_minus_HeI_feedback(line, znow, XHI1s, HeIA);
        
    //====================================================================
    // correction for case Pd~1 and pd != 1
    //====================================================================
    double Pd=p_ij(pd*tauS), psc=1.0-pd;
    double fc =psc*Pd/(1.0-psc*Pd);
    if(pd*tauS<=1.0e-6) fc=psc/pd*(1.0-tauS*pd/2.0)/(1.0+psc*tauS/2.0);
        
    HeI_photons[line].z[iznow]=znow;
    //====================================================================
    // these are the distortions at redshift z!
    //====================================================================
    if(pd*tauS<=1.0e-6)
    { 
        HeI_photons[line].Dn[iznow]=(DnL*pd+fc*pd*(DnL-Dn))*tauS*(1.0-tauS*pd/2.0); 
    }
    else
    {
        HeI_photons[line].Dn[iznow]=DnL+fc*(DnL-Dn);
        //================================================================
        // for line with taud<~1 one will have a factor of 1-e^-taud < 1
        //================================================================
        HeI_photons[line].Dn[iznow]*=1.0-exp(-tauS*pd);
    }
    
    if(mess!="" && He_feedback_message) cout << " update_current_feedback_functions " 
         << mess << " : " << line << " i= " << " " <<  HeI_Lines[line].upperLevel.i << " z= " 
         << znow << " (" << HeI_Lines[line].upperLevel.n << ", " << HeI_Lines[line].upperLevel.l 
         << ", " << HeI_Lines[line].upperLevel.s << ", " << HeI_Lines[line].upperLevel.j 
         << ") tauS= " << tauS << " pd= " << pd << " " << HeI_photons[line].Dn[iznow] 
         << "\n\t\t DnL*(1-e^-taud)= " << DnL*(1.0-exp(-tauS*pd)) << " fc= " << fc
         << " ( DI(z=0)= " << const_h*nu21/znow*2.0*pow(nu21/const_cl/znow, 2)*HeI_photons[line].Dn[iznow] 
         << ", nu(z=0)= " << nu21/znow/1.0e+9 <<  " GHz) " 
         << Dn << " " << Dn*exp(-pd*tauS) 
         << " sum= " << HeI_photons[line].Dn[iznow]+Dn*exp(-pd*tauS) << endl;
    
    //====================================================================
    // part from line that is passing through in addition
    //====================================================================
    HeI_photons[line].Dn[iznow]+=Dn*exp(-pd*tauS);
    
    return;
}

//==================================================================================================
void save_HeI_distortion(Gas_of_HeI_Atoms &HeIA, double XHI1s, int iznow, double znow)
{
    //====================================================================
    // save information for HeIS-transition
    //====================================================================
    HeI_photons_iznow=iznow;
    
    //====================================================================
    // update the other way around
    //====================================================================
    if(He_feedback_message) cout << endl;

    for(int line=(int)HeI_Lines.size()-1; line>=1; line--) 
        update_current_feedback_functions(line, HeIA, XHI1s, iznow, znow, "save_HeI_distortion");
        
    return;
}

//==================================================================================================
//==================================================================================================
