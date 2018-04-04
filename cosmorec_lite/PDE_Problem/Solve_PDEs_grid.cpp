//====================================================================================================================
// setup for isolated resonances; if the upper/lower boundary comes too close a linear grid will be used
//====================================================================================================================
int init_PDE_xarr_HI(double *xarr, double xmin, double xmax, double xres, double xD, int mess=0)
{
    xD*=xres;

    //===================================================================================
    // parameters for grid setup
    //===================================================================================
    double core_width=10.0;    // core width in Doppler units
    double enhance_fac_up=2.4; // enhancement of points above resonance
    
    //===================================================================================
    // number of grid-points in core
    //===================================================================================
    double Doppler_core=(min(xmax, xres+core_width*xD)-max(xmin, xres-core_width*xD))/xD; 
    int ncore=core_HI*Doppler_core;
    
    //===================================================================================
    // number of Doppler width above and below the core region. If the lines become very 
    // close to eachother then a linear grid is used in the wing region. For this a 
    // minimal density of points per Doppler width is set.
    //===================================================================================
    double distance_treshold=40.0;
    double Doppler_up=(xmax-xres)/xD-core_width;
    double Doppler_low=(xres-xmin)/xD-core_width;

    int diluted_core=core_HI/2; 
    int init_up=1, init_low=1;
    double min_np_up=10, min_np_low=10;
    
    if(Doppler_up <distance_treshold){ init_up=0;  min_np_up =Doppler_up *diluted_core; }
    if(Doppler_low<distance_treshold){ init_low=0; min_np_low=Doppler_low*diluted_core; }
    
    if(mess==1) cout << " xres= " << xres << endl;
    
    //===================================================================================
    double dec_up=log10(xmax-xres)-log10(core_width*xD);
    double dec_low=log10(xres-xmin)-log10(core_width*xD);
    
    //===================================================================================   
    int nup=(int)max(min_np_up, enhance_fac_up*dec_HI*dec_up);
    int nlow=(int)max(min_np_low, dec_HI*dec_low);
    
    //===================================================================================   
    if(Doppler_up<=0.0) nup=0;
    if(Doppler_low<=0.0) nlow=0;
    
    //===================================================================================   
    int nused=nup+nlow+ncore;

    //===================================================================================
    vector<double> Dxarr(nused);

    if(nup>0)
    {
        init_xarr(core_width*xD, xmax-xres, &xarr[nlow+ncore], nup, init_up, mess);
        for(int k=0; k<nup; k++) xarr[nlow+ncore+k]+=xres;
    }

    init_xarr(max(xmin, xres-core_width*xD), min(xmax, xres+core_width*xD), &xarr[nlow], ncore+1, 0, mess);

    if(nlow>0)
    {
        init_xarr(core_width*xD, xres-xmin, &Dxarr[0], nlow+1, init_low, mess);
        for(int k=0; k<nlow+1; k++) xarr[k]=xres-Dxarr[nlow-k];
    }
    
    for(int k=1; k<nused; k++) if(xarr[k-1]>xarr[k])
    {
        cerr << " init_PDE_xarr_HI:: grid non-monotonic " << k-1 << " " << xarr[k-1] << endl;
        cerr << "                                        " << k << " " << xarr[k] << endl; exit(0);
    }

    return nused;
}

//====================================================================================================================
// arrangement of grid
//====================================================================================================================
void init_PDE_xarr_cores_HI(vector<double> &xarr, double xmin, double xmax, vector<double> &resonances)
{
    xarr.resize(100000); // just make the array (very) large for the moment; unused points will be deleted
    
    if(show_messages>=1) cout << "\n init_PDE_xarr_cores_HI :: setting up grid. " << endl;
    
    int mess=(show_messages>=2 ? 1 : 0);
    int nres=resonances.size();
    int np_2s=npts_2s1s; 
    double xcrit=0.90;
    double xD_f_xarr = 2.0e-5;    

    //=================================================================================
    // introduce small region around upper boundary that is treated with linear grid
    //=================================================================================
    double xmax_all=xmax;
    xmax=max(xmax*0.98, 0.75*xmax+0.25*resonances[nres-1]);
    
    vector<double> dx(nres+1);
    dx[0]=xmax-xmin;
    
    //=================================================================================
    // define regions around resonances
    //=================================================================================
    if(nres>1) dx[0]=(resonances[0]+resonances[1])*0.5-xmin;
    for(int i=1; i<nres-1; i++) dx[i]=(resonances[i+1]-resonances[i-1])*0.5;
    if(nres>1) dx[nres-1]=xmax-(resonances[nres-2]+resonances[nres-1])*0.5;
    
    //=================================================================================
    // below Ly-n resonances (2s-1s part)
    //=================================================================================
    double xl=xmin;
    int istart=0;

    if(xmin<=xcrit)
    {
        double xcrit_2=0.05;        
        int np_2s_log=0;
        if(xmin<=xcrit_2)
        {
            int dens_log=30;
            double dlog=log10(xcrit_2/xmin);
            np_2s_log=dlog*dens_log;
            
            //=============================================================================
            // normal log-scale for x<xcrit_2
            //=============================================================================
            init_xarr(xl, xcrit_2, &xarr[istart], np_2s_log, 1, mess);    
            istart+=np_2s_log-1;
            xl=xcrit_2;
        }
            
        //=============================================================================
        // linear scale between xcrit_2 and xcrit
        //=============================================================================
        init_xarr(xl, xcrit, &xarr[istart], np_2s-np_2s_log, 0, mess);    
        istart+=np_2s-np_2s_log-1;
        xl=xcrit;
        
        dx[0]+=xmin-xl;
    }
    
    if(nres>1) 
    {
        //=============================================================================
        // Ly-a
        //=============================================================================
        int nused=init_PDE_xarr_HI(&xarr[istart], xl, xl+dx[0], resonances[0], xD_f_xarr, mess);
        istart+=nused-1;
        
        //=============================================================================
        // Ly-n
        //=============================================================================
        for(int i=1; i<nres; i++) 
        {
            xl+=dx[i-1];
            nused=init_PDE_xarr_HI(&xarr[istart], xl, xl+dx[i], resonances[i], xD_f_xarr, mess);
            istart+=nused-1;
        }
    }
    else
    {
        //=============================================================================
        // just one resonance present
        //=============================================================================
        int nused=init_PDE_xarr_HI(&xarr[istart], xl, xl+dx[0], resonances[0], xD_f_xarr, mess);
        istart+=nused-1;
    }
    
    //=================================================================================
    // close off with linear grid
    //=================================================================================
    int nused=100;
    if(nres==0) nused=npts_2s1s;
    init_xarr(xmax, xmax_all, &xarr[istart], nused, 0, mess);    
    istart+=nused-1;
    
    for(int k=1; k<istart; k++) if(xarr[k-1]>xarr[k])
    { cerr << " init_PDE_xarr_cores_HI::grid non-monotonic " << k << " " << xarr[k] << endl; exit(0); }
    
    if(show_messages>=1) cout << " init_PDE_xarr_cores_HI:: done. " << endl;

    //=================================================================================
    // deleted unused points.
    //=================================================================================
    xarr.erase(xarr.begin()+istart+1, xarr.end());

    return;
}   


