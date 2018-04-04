//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// This module allows to add the recombination corrections as described by 
// Chluba & Thomas 2010.
//========================================================================================

#ifndef REC_CORRS_CT_H
#define REC_CORRS_CT_H

//========================================================================================
// to load the correction factor
//========================================================================================
void read_recombination_correction(int mflag=1);

//========================================================================================
// correction factor 1+DNe/Ne as explained by Rubino-Martin et al. 2010.
//========================================================================================
double Recombination_correction_factor(double z);

#endif
