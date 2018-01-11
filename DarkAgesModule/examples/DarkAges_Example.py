# -*- coding: utf-8 -*-

import sys
import os
import numpy as np

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.insert(1,root_path)

from DarkAges import get_redshift, get_logEnergies, transfer_functions, DarkAgesError, DarkOptions, channel_dict
from DarkAges.common import sample_spectrum, finalize
from DarkAges.model import annihilating_model as model

#######################
#######################
#######################

specfile = sys.argv[1]

# spectrum of 'dN/dE'-type
spec_data = np.genfromtxt(specfile, unpack=True, usecols=(0,1,2,3,4), skip_header=1, dtype=np.float64)
masses = np.unique(spec_data[0])
try:
	assert len(masses) == 1
	mass = masses[0]
except AssertionError:
	raise DarkAgesError('It seems that the file {:s} contains spectra for several masses. (Found {:d} different masses)'.format(specfile,len(masses)))

logEnergies = get_logEnergies()
redshift = get_redshift()

spec_el, spec_ph, spec_oth = sample_spectrum(spec_data[2], spec_data[3], spec_data[4], spec_data[1], mass, logEnergies, **DarkOptions)
example_model = model(spec_el, spec_ph, spec_oth, 1e9*mass, logEnergies,redshift)

f_function = np.zeros( shape=(5,len(redshift)), dtype=np.float64 )
for channel in channel_dict:
	idx = channel_dict[channel]
	f_function[idx,:] = example_model.calc_f(transfer_functions[idx])[-1]

### Finalize (Print the output to be read by CLASS) ###

finalize(redshift,
         f_function[channel_dict['Heat']],
         f_function[channel_dict['Ly-A']],
         f_function[channel_dict['H-Ion']],
         f_function[channel_dict['He-Ion']],
         f_function[channel_dict['LowE']],
         **DarkOptions)
