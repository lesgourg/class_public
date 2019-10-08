from __future__ import absolute_import, division, print_function

import numpy as np
import os
import sys
import dill
from scipy.special import erf

from .__init__ import DarkAgesError as err
data_dir = os.path.join( os.path.dirname(os.path.realpath( __file__ )), 'data' )

def boost_factor_halos(redshift,zh,fh):
	ret = 1 + fh*erf(redshift/(1+zh))/redshift**3
	return ret

def secondaries_from_cirelli(logEnergies,mass,primary, **DarkOptions):
	from .common import sample_spectrum
	cirelli_dir = os.path.join(data_dir, 'cirelli')
	dumpername = 'cirelli_spectrum_of_{:s}.obj'.format(primary)

	injection_history = DarkOptions.get("injection_history","annihilation")
	if "decay" in injection_history:
		equivalent_mass = mass/2.
	else:
		equivalent_mass = mass
	if equivalent_mass < 5 or equivalent_mass > 1e5:
		raise err('The spectra of Cirelli are only given in the range [5 GeV, 1e2 TeV] assuming DM annihilation. The equivalent mass for the given injection_history ({:.2g} GeV) is not in that range.'.format(equivalent_mass))

	if not hasattr(logEnergies,'__len__'):
		logEnergies = np.asarray([logEnergies])
	else:
		logEnergies = np.asarray(logEnergies)

	if not os.path.isfile( os.path.join(cirelli_dir, dumpername)):
		sys.path.insert(1,cirelli_dir)
		from spectrum_from_cirelli import get_cirelli_spectra
		masses, log10X, dNdLog10X_el, dNdLog10X_ph, dNdLog10X_oth = get_cirelli_spectra(primary)
		total_dNdLog10X = np.asarray([dNdLog10X_el, dNdLog10X_ph, dNdLog10X_oth])
		from .interpolator import NDlogInterpolator
		interpolator = NDlogInterpolator(masses, np.rollaxis(total_dNdLog10X,1), exponent = 0, scale = 'log-log')
		dump_dict = {'dNdLog10X_interpolator':interpolator, 'log10X':log10X}
		with open(os.path.join(cirelli_dir, dumpername),'wb') as dump_file:
			dill.dump(dump_dict, dump_file)
	else:
		with open(os.path.join(cirelli_dir, dumpername),'rb') as dump_file:
			dump_dict = dill.load(dump_file)
			interpolator = dump_dict.get('dNdLog10X_interpolator')
			log10X = dump_dict.get('log10X')
	del dump_dict
	temp_log10E = log10X + np.log10(equivalent_mass)*np.ones_like(log10X)
	temp_el, temp_ph, temp_oth = interpolator.__call__(equivalent_mass) / (10**temp_log10E * np.log(10))[None,:]
	ret_spectra = np.empty(shape=(3,len(logEnergies)))
	ret_spectra = sample_spectrum(temp_el, temp_ph, temp_oth, temp_log10E, mass, logEnergies, **DarkOptions)
	return ret_spectra

def secondaries_from_simple_decay(E_secondary, E_primary, primary):
	if primary not in ['muon','pi0','piCh']:
		raise err('The "simple" decay spectrum you asked for (species: {:s}) is not (yet) known.'.format(primary))

	if not hasattr(E_secondary,'__len__'):
		E_secondary = np.asarray([E_secondary])
	else:
		E_secondary = np.asarray(E_secondary)

	decay_dir  = os.path.join(data_dir, 'simple_decay_spectra')
	dumpername = 'simple_decay_spectrum_of_{:s}.obj'.format(primary)
	original_data = '{:s}_normed.dat'.format(primary)

	if not os.path.isfile( os.path.join(decay_dir, dumpername)):
		data = np.genfromtxt( os.path.join(decay_dir, original_data), unpack = True, usecols=(0,1,2,3))
		from .interpolator import NDlogInterpolator
		spec_interpolator = NDlogInterpolator(data[0,:], data[1:,:].T, exponent = 1, scale = 'lin-log')
		dump_dict = {'spec_interpolator':spec_interpolator}
		with open(os.path.join(decay_dir, dumpername),'wb') as dump_file:
			dill.dump(dump_dict, dump_file)
	else:
		with open(os.path.join(decay_dir, dumpername),'rb') as dump_file:
			dump_dict = dill.load(dump_file)
			spec_interpolator = dump_dict.get('spec_interpolator')

	x = E_secondary / E_primary
	out = spec_interpolator.__call__(x)
	out /= (np.log(10)*E_secondary)[:,None]
	return out

def luminosity_accreting_bh(Energy,recipe,PBH_mass):
	if not hasattr(Energy,'__len__'):
		Energy = np.asarray([Energy])
	if recipe=='spherical_accretion':
		a = 0.5
		Ts = 0.4*511e3
		Emin = 1
		Emax = Ts
		out = np.zeros_like(Energy)
		Emin_mask = Energy > Emin
		# Emax_mask = Ts > Energy
		out[Emin_mask] = Energy[Emin_mask]**(-a)*np.exp(-Energy[Emin_mask]/Ts)
		out[~Emin_mask] = 0.
		# out[~Emax_mask] = 0.

	elif recipe=='disk_accretion':
		a = -2.5+np.log10(PBH_mass)/3.
		Emin = (10/PBH_mass)**0.5
		# print a, Emin
		Ts = 0.4*511e3
		out = np.zeros_like(Energy)
		Emin_mask = Energy > Emin
		out[Emin_mask] = Energy[Emin_mask]**(-a)*np.exp(-Energy[Emin_mask]/Ts)
		out[~Emin_mask] = 0.
		Emax_mask = Ts > Energy
		out[~Emax_mask] = 0.
	else:
		from .__init__ import DarkAgesError as err
		raise err('I cannot understand the recipe "{0}"'.format(recipe))
	# print out, Emax_mask
	return out/Energy #We will remultiply by Energy later in the code
