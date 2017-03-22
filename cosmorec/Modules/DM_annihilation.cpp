//========================================================================================
// Author: Jens Chluba (Nov. 19th 2010)
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// This simple module allows to add the ionizing effect of DM-annihilation. The important 
// equations as implemented here were given in Chluba 2010. Original works are 
// Chen & Kamionkowski 2004 and Padmanabhan & Finkbeiner 2005.
// In this module additional excitations are not included.
//========================================================================================

//========================================================================================
// branching into ionizations according to Chen & Kamionkowski 2004
// (see also Chluba 2010)
//========================================================================================
double branching_ratio_ions_Chen(double x){ return (1.0-x)/3.0; }

double branching_ratio_heat_Chen(double xp, double ZHeII, double fHe)
{ return (1.0 + 2.0*xp + fHe*(1.0 + 2.0*ZHeII))/3.0/(1.0+fHe); }


//========================================================================================
// branching into ionizations according to Shull & van Steenberg 1985 
// (see also Chluba 2010)
//========================================================================================
double branching_ratio_ions_SvSt(double x)
{ return 0.4*pow(1.0-pow(x, 0.4), 7.0/4.0); }


//========================================================================================
// evalutation of ODE terms
//========================================================================================
void evaluate_DM_annihilation_terms(double z, double fHe, double Tg,
                                    double XHII, double XHeII, double fDM, 
                                    double &dXHI_1s_dt, double &dXHeI_1s_dt, 
                                    double &drho_dt)
{
    //====================================================================================
    // fDM [eV/s] gives annihilation efficiency; typical value fDM=2.0e-24 eV/s
    // (see Chluba 2010 for definitions)
    // In this module additional excitations are not included at the moment
    //====================================================================================
    
    double Fcal=pow(1.0+z, 3);
    double dE_dt=fDM*Fcal;  // here the factor 1/NH was canceled with definition of dE_dt
    //
    double bHe_ion=branching_ratio_ions_Chen(XHeII/fHe);
    double bHI_ion=branching_ratio_ions_Chen(XHII);
    double b_heat =branching_ratio_heat_Chen(XHII, XHeII/fHe, fHe);
    //
    double fheat=2.0/3.0*const_e/const_kB/(1.0+fHe+XHII+XHeII);
    
    dXHI_1s_dt+=-bHI_ion/13.6/(1.0+fHe)*dE_dt;
    if(XHeII>0.0) dXHeI_1s_dt+=-bHe_ion/24.6*fHe/(1.0+fHe)*dE_dt;
    drho_dt+=fheat*b_heat*dE_dt/Tg;
    
    return;
}
