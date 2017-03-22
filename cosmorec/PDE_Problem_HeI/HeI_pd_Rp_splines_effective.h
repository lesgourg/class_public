//===========================================================================================================
// Authors Jeffrey Fung & Jens Chluba Feb/March 2010
//
// comment: HeI populations have to be pre-loaded; Also the effective rates have to be initialized
//          Also note that here calc_Xi and all other splines are using the resolved level index 
//          instead of the HeI-level index!
//===========================================================================================================
// 30.07.2014: Xi_Data is passed on read only
// 22.07.2014: restart of spline setup if number of redshift points changes
// 22.07.2014: changed setup routine for splines (smaller buffer region)

#ifndef HEI_PD_RP_SPLINES_EFFECTIVE_H
#define HEI_PD_RP_SPLINES_EFFECTIVE_H

#include <string>
#include "HeI_Atom.h"
#include "Cosmos.h"

//===========================================================================================================
// For precalulated death probabilities and emission rates
//===========================================================================================================
void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_HeI_Atoms &HeA);

void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_HeI_Atoms &HeA,
                                            const vector<vector<double> > &X_Data);

void clear_HeI_pd_Rp_effective_memory();
void Set_HeI_pd_Rp_splines_effective_verbosity(int v);

//===========================================================================================================
// pd-functions (the index ires is the corresponding resolved level index)
//===========================================================================================================
double calc_HeI_pd_i_splines_effective(double z, int ires);

//===========================================================================================================
// Dnem-functions (the index ires is the corresponding resolved level index)
//===========================================================================================================
double calc_HeI_Dnem_i_splines_effective(double z, int ires);

#endif

//===========================================================================================================
//===========================================================================================================
