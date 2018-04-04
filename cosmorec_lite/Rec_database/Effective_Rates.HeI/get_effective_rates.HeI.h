//======================================================================
// get_effective_rates.HeI.cpp
//
// Created by Jens Chluba (Oct 2010).
// Copyright 2010 CITA (U of T). All rights reserved.
//
//======================================================================

#ifndef GET_EFFECTIVE_RATES_HEI_H
#define GET_EFFECTIVE_RATES_HEI_H

#include <string>
#include <vector>

#include "HeI_Atom.h"

void read_parameter_file(string fname, Gas_of_HeI_Atoms &HeIA, vector<int> &res_Level_indices);

void load_rates(string path, int nShell_max, Gas_of_HeI_Atoms &HIA);

void get_rates_HeI(double Tg, vector<double> &Ai, vector<double> &Bi, 
                   vector<vector<double> > &RijVec);

void get_rates_all_HeI(double Tg, vector<double> &Ai, vector<double> &Bi, 
                       vector<double> &Bitot, vector<vector<double> > &RijVec);

void get_Bitot_HeI(double Tg, vector<double> &Bitot);
void get_Bitot_HeI(double Tg, double &Bitot, int ires);

int get_HeI_index(int index_of_res_state);
int get_number_of_resolved_levels_HeI();
int get_res_index_HeI(int index_of_HeI_state);

#endif
