//===========================================================================================================
// Author Jens Chluba May 2009
// comment: HI populations have to be pre-loaded; Also the effective rates have to be initialized
//===========================================================================================================
// 30.07.2014: Xi_Data is passed on read only
// 22.07.2014: restart of spline setup if number of redshift points changes
// 22.07.2014: changed setup routine for splines (smaller buffer region)

#ifndef HI_PD_RP_SPLINES_EFFECTIVE_H
#define HI_PD_RP_SPLINES_EFFECTIVE_H

#include <string>
#include "Atom.h"
#include "Cosmos.h"

//===========================================================================================================
// For precalulated death probabilities and emission rates
//===========================================================================================================
void set_up_splines_for_HI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_Atoms &HA);

void set_up_splines_for_HI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_Atoms &HA, 
                                           const vector<vector<double> > &X_Data);

void clear_HI_pd_Rp_effective_memory();
void Set_HI_pd_Rp_splines_effective_verbosity(int v);

//===========================================================================================================
// pd-functions
//===========================================================================================================
double calc_HI_pd_nl_splines_effective(double z, int n, int l);
double calc_HI_pd_i_splines_effective(double z, int i);

//===========================================================================================================
// Dnem-functions
//===========================================================================================================
double calc_HI_Dnem_nl_splines_effective(double z, int n, int l);
double calc_HI_Dnem_i_splines_effective(double z, int i);

#endif

//===========================================================================================================
//===========================================================================================================
