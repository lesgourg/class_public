from __future__ import absolute_import, division, print_function

import numpy as np
import os

column_dict = {
	'electron': 4,
	'muon': 7,
	'tau':	10,
	'quark':11,
	'charm':12,
	'bottom':13,
	'top':14,
	'wboson':17,
	'zboson':20,
	'gluon':21,
	'photon':22,
	'higgs':23
}

data_dir = os.path.dirname(os.path.realpath( __file__ ))

def get_cirelli_spectra(key):
	EW_cols = (0, 1, column_dict.get(key))

	data_elec_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_positrons.dat'), unpack=True, usecols=EW_cols, skip_header=1)
	data_phot_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_gammas.dat'), unpack=True, usecols=EW_cols, skip_header=1)
	data_nu_e_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_neutrinos_e.dat'), unpack=True, usecols=EW_cols, skip_header=1)
	data_nu_m_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_neutrinos_mu.dat'), unpack=True, usecols=EW_cols, skip_header=1)
	data_nu_t_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_neutrinos_tau.dat'), unpack=True, usecols=EW_cols, skip_header=1)
	data_prot_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_antiprotons.dat'), unpack=True, usecols=EW_cols, skip_header=1)
	data_deut_EW = np.genfromtxt(os.path.join(data_dir,'AtProduction_antideuterons.dat'), unpack=True, usecols=EW_cols, skip_header=1)

	masses = np.unique(data_elec_EW[0,:])
	log10X = np.unique(data_elec_EW[1,:])
	dim1 = len(masses)
	dim2 = len(log10X)
	dNdlog10X_el = 2*data_elec_EW[2,:].reshape(dim1,dim2)
	dNdlog10X_ph = data_phot_EW[2,:].reshape(dim1,dim2)
	dNdlog10X_oth = 2*(data_nu_e_EW[2,:] + data_nu_m_EW[2,:] + data_nu_t_EW[2,:] + data_prot_EW[2,:] + data_deut_EW[2,:]).reshape(dim1,dim2)

	return masses, log10X, dNdlog10X_el, dNdlog10X_ph, dNdlog10X_oth 

