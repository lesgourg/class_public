import itertools
import dill
import numpy as np
import sys
import os 
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges.common import channel_dict, get_index
from DarkAges import logEnergies, redshift, transfer_functions
from DarkAges.recipes import load_from_spectrum
from DarkAges.model import model
from DarkAges.interpolator import logInterpolator

model_dir = os.path.split(os.path.realpath(__file__))[0]
model_name =  model_dir.split('/')[-1]

masses = np.array([50,55,60,65,70,75,80,85,90,95,100], dtype=np.int32)
primaries = np.array(['muon_StandAlone','bottom'], dtype=np.dtype('a32'))
f_grid = np.empty(shape=(len(masses),len(primaries),5,len(redshift)), dtype = np.float64)

for mass, primary in itertools.product(*(masses,primaries)):
	idx_mass = get_index(masses, mass)
	idx_prim = get_index(primaries, primary)
	fname = os.path.join(model_dir, 'data/{:s}_{:d}_spectrum.dat'.format(primary,mass))
	tmp_model = load_from_spectrum(fname, logEnergies, mass, np.inf)
	for channel in channel_dict:
		idx_chan = channel_dict[channel]
		f_grid[idx_mass,idx_prim,idx_chan,:] = tmp_model.calc_f(transfer_functions[idx_chan])[-1]
		
interpolated_f = np.empty(shape=(len(primaries),5,len(redshift)), dtype= logInterpolator)
for primary, idx_chan, z in itertools.product(*(primaries,channel_dict.values(),redshift)):
	idx_prim = get_index(primaries, primary)
	idx_z = get_index(redshift, z)
	interpolated_f[idx_prim, idx_chan, idx_z] = logInterpolator(masses, f_grid[:,idx_prim,idx_chan,idx_z], 0)	

dump_dict = {'f-function': interpolated_f}

with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'wb') as dump_file:
	dill.dump(dump_dict, dump_file)
	
		
		


