//============================================================================
// Author: Jens Chluba
// July 2007
// Purpose: compute the correction to the escape probability in the helium 
// lines caused by the absorption from hydrogen
//============================================================================
#ifndef PESC_HI_H
#define PESC_HI_H

#include "Photoionization_cross_section.h"
#include "Voigtprofiles.h"

//===================================================================================
// Simple approximations for DPesc
//===================================================================================
// Inner integral performed analtically & symmetry is used
double DPesc_appr_I_sym(double z, double eta_S, double eta_c, Voigtprofile_Dawson &P, 
                        Photoionization_cross_section_SH &L, double abs=0.0);

// Inner integral performed analtically & symmetry and log coordinates are used
double DPesc_appr_I_sym_lg(double z, double eta_S, double eta_c, Voigtprofile_Dawson &P, 
                           Photoionization_cross_section_SH &L, double abs=0.0);

#endif
