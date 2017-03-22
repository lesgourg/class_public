//===========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//===========================================================================================

//===========================================================================================
// This simple module allows to add the ionizing effect of DM-annihilation. The important 
// equations as implemented here were given in Chluba 2010. Original works are 
// Chen & Kamionkowski 2004 and Padmanabhan & Finkbeiner 2005.
// In this module additional excitations are not included at the moment.
//===========================================================================================

#ifndef DM_ANNIHILATION_H
#define DM_ANNIHILATION_H

void evaluate_DM_annihilation_terms(double z, double Hz, double fHe, double xp, double xHep, 
                                    double CHe, double CH, double fDM, double fvec[]);

#endif
