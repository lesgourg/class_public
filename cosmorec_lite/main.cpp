//==================================================================================================
// Author: Jens Chluba 
// last modification: Dec 2010
// purpose: to access the main CosmoRec code from console
//==================================================================================================

#include "routines.h"
#include "CosmoRec.h"

#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

//==================================================================================================
//
// to check CosmoRec over a range of different cosmologies
//
//==================================================================================================
int main_random()
{ 
    int nrep = 100000;    
    int runmode=0;
    // 0: DM annihilation; 1: run settings; 2: outputs; 3: A2s1s; 4: dummy
    double runpars[5]={0.0, 0, 2, 0.0, 0.0};
    
    //==============================================================================================
    double ombh2_l = 0.01;
    double ombh2_u = 0.04;
    
    double omch2_l = 0.05;
    double omch2_u = 0.22;
    
    double h0_l = 0.6;
    double h0_u = 0.8;

    double T0_l = 2.6;
    double T0_u = 2.8;    

    double yhe_l = 0.1;
    double yhe_u = 0.6;    

    double Nnu_l = 1.0;
    double Nnu_u = 6.0;    

    double A2s1s_l = 4.0;
    double A2s1s_u = 15.0;
    //==============================================================================================
    
    double h0, omegab, omegac, omegak=0.0;
    double tcmb, yhe, Nnu, A2s1s;
    
    int nz=10000;
    double *z_arr=new double[nz];
    
    init_xarr(0.01, 1.0e+4, z_arr, nz, 0, 0);
    
    srand48(time(0));
    
    double *Xe_arr=new double[nz];
    double *Te_arr=new double[nz];
    
    cout.precision(10); 
    clock_t seconds_start=clock(), seconds_end;

    for(int i = 0; i < nrep; i++) 
    {
        cout << " **** Iteration: " << i << " ****" << endl;

        h0 = h0_l + (h0_u - h0_l) * drand48();
        omegab = (ombh2_l + (ombh2_u - ombh2_l) * drand48()) / (h0*h0);
        omegac = (omch2_l + (omch2_u - omch2_l) * drand48()) / (h0*h0);
        tcmb = T0_l + (T0_u - T0_l) * drand48();
        yhe = yhe_l + (yhe_u - yhe_l) * drand48();
        Nnu = Nnu_l + (Nnu_u - Nnu_l) * drand48();
        A2s1s = A2s1s_l + (A2s1s_u - A2s1s_l) * drand48();
      
        cout << " omegab= " << omegab << endl;
        cout << " omegam= " << omegac+omegab << endl;
        cout << " h0=  " << h0 << endl;
        cout << " tcmb=  " << tcmb << endl;
        cout << " yhe=  " << yhe << endl;
        cout << " Nnu=  " << Nnu << endl;
        cout << " A2s1s=  " << A2s1s << endl;
        
        runpars[3]=A2s1s;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu,
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, 0);
        
        seconds_end=clock();
        double seconds_Drun=(1.0*seconds_end-seconds_start)/CLOCKS_PER_SEC;
        cout << " The run took " << seconds_Drun << " sec. " << endl;    
        seconds_start=seconds_end;

        cout << " **** Done. ****" << endl << endl;    
    }

    //===============================================================
    // clean up 
    //===============================================================
    delete [] z_arr;
    delete [] Xe_arr;
    delete [] Te_arr;
    
    cleanup_CosmoRec();
    
    return 0;
}

//==================================================================================================
//
// example how to call CosmoRec in batch-mode
//
//==================================================================================================
int main_batch_mode()
{ 
    int runmode=0, label=1;
    double runpars[5]={0.0, 0, 0, 0.0, 0.0};
    
    double omegac=0.216;
    double omegab=0.044;
    double Nnu=3.04;
    double omegak=0.0;
    double h0=0.71;
    double tcmb=2.725;
    double yhe=0.24;
    
    int nz=10000;
    double *z_arr=new double[nz];
    double *Xe_arr=new double[nz];
    double *Te_arr=new double[nz];
    
    init_xarr(0.01, 1.0e+4, z_arr, nz, 0, 0);
    
    //===============================================================
    // runs
    //===============================================================
    int count=0;
    for(int k=0; k<2; k++)
    {
        tcmb=2.725;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.8; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.65; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.75; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.7; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.764; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.766; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        tcmb=2.743; label++;
        CosmoRec(runmode, runpars, omegac, omegab, omegak, Nnu, 
                 h0, tcmb, yhe, nz, z_arr, Xe_arr, Te_arr, label);
        
        count+=8;
    }
    
    cout << count << endl;
    
    //===============================================================
    // clean up 
    //===============================================================
    delete [] z_arr;
    delete [] Xe_arr;
    delete [] Te_arr;
    
    cleanup_CosmoRec();
    
    return 0;
}

//==================================================================================================
// main 
//==================================================================================================
int main(int narg, char *args[])
{ 
    if(narg>=2) return CosmoRec(narg, args);
//    return main_batch_mode();
    return main_random();
}
