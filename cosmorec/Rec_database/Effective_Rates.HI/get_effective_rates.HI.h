//======================================================================
// get_effective_rates.HI.h
//
//
// Created by Rajat Thomas and Jens Chluba on 10-07-20.
// Copyright 2010 CITA ( U of T). All rights reserved.
//
//======================================================================

#ifndef GET_EFFECTIVE_RATES_H
#define GET_EFFECTIVE_RATES_H

#include <string>
#include <vector>

#include "Atom.h"

void load_rates(string path, int nShell_max, Gas_of_Atoms &HIA);

void get_rates(double Tg, double Te, vector<double> &Ai, vector<double> &Bi, 
               vector<vector<double> > &RijVec);

void get_rates_all(double Tg, double Te, vector<double> &Ai, vector<double> &Bi, 
                   vector<double> &Bitot, vector<vector<double> > &RijVec);

void get_Bitot(double Tg, vector<double> &Bitot);
void get_Bitot(double Tg, double &Bitot, int ires);

int get_HI_index(int index_of_res_state);
int get_number_of_resolved_levels();
int get_res_index(int index_of_HI_state);

#endif
