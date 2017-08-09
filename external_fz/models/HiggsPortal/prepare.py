import itertools
import dill
import numpy as np
import sys
import os 
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges.common import channel_dict, get_index
from DarkAges import logEnergies, redshift, transfer_functions, options
from DarkAges.recipes import load_from_spectrum
from DarkAges.model import model
from DarkAges.interpolator import logInterpolator, logLinearInterpolator

model_dir = os.path.split(os.path.realpath(__file__))[0]
model_name =  model_dir.split('/')[-1]

masses = np.genfromtxt(os.path.join(model_dir, 'data/Masses.dat'), unpack=True, usecols=(0), dtype=np.float64)
masses = masses[masses <= 5e3]
frac_read = np.genfromtxt(os.path.join(model_dir,'data/Fractions.dat'), unpack=True, usecols=(0,4,5,6,7,8,9,10,11,12), dtype=np.float64)
primaries = np.array(['top', 'higgs', 'zzstar', 'wwstar', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a32'))
#primaries = np.array(['top', 'higgs'], dtype=np.dtype('a32'))

fractions_interp_obj = np.empty(shape=(len(primaries)), dtype=logLinearInterpolator)
for idx in xrange(len(primaries)):
	fractions_interp_obj[idx] = logLinearInterpolator(frac_read[0,:], frac_read[idx+1,:], 0)

total = len(masses) * len(primaries)
counter = 0. 

f_grid = np.empty(shape=(len(masses),len(primaries),5,len(redshift)), dtype = np.float64)
for mass, primary in itertools.product(*(masses,primaries)):
	counter += 1.
	print (counter / total) * 100.
	idx_mass = get_index(masses, mass)
	idx_prim = get_index(primaries, primary)
	fname = os.path.join(model_dir, 'data/{:s}_EW.dat'.format(primary))
	tmp_model = load_from_spectrum(fname, logEnergies, mass, np.inf)
	for channel in channel_dict:
		idx_chan = channel_dict[channel]
		f_grid[idx_mass,idx_prim,idx_chan,:] = tmp_model.calc_f(transfer_functions[idx_chan])[-1]
		
interpolated_f = np.empty(shape=(len(primaries),5,len(redshift)), dtype= logInterpolator)
for primary, idx_chan, z in itertools.product(*(primaries,channel_dict.values(),redshift)):
	idx_prim = get_index(primaries, primary)
	idx_z = get_index(redshift, z)
	interpolated_f[idx_prim, idx_chan, idx_z] = logInterpolator(masses, f_grid[:,idx_prim,idx_chan,idx_z], 0)	

dump_dict = {'fractions':fractions_interp_obj, 'f-functions': interpolated_f}

with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'wb') as dump_file:
	dill.dump(dump_dict, dump_file)
	
		
		


