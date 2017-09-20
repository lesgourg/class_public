import dill
import numpy as np
import sys
import os
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges.common import channel_dict
from DarkAges import logEnergies
from DarkAges.recipes import load_from_spectrum
from DarkAges.interpolator import logInterpolator, logLinearInterpolator

model_dir = os.path.split(os.path.realpath(__file__))[0]
model_name =  model_dir.split('/')[-1]

def prepare():
	dump_dict = dict()

	frac_read = np.genfromtxt(os.path.join(model_dir,'data/Fractions.dat'), unpack=True, usecols=(0,4,5,6,7,8,9,10,11,12), dtype=np.float64)
	primaries = np.array(['top', 'higgs', 'zzstar', 'wwstar', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a32'))
	#primaries = np.array(['top', 'higgs'], dtype=np.dtype('a32'))

	dump_dict.update({'channels':primaries})

	for idx_prim, primary in enumerate(primaries):
		dump_dict_key = 'frac_interp_{:s}'.format(primary)
		dump_dict_value = logLinearInterpolator(frac_read[0,:], frac_read[idx_prim+1,:], 0,logspace_wanted=False)
		dump_dict.update({dump_dict_key:dump_dict_value})

		fname = os.path.join(model_dir, 'data/{:s}_EW.dat'.format(primary))
		dump_dict_key = 'spec_interp_{:s}'.format(primary)
		dump_dict_value = load_from_spectrum(fname, logEnergies)
		dump_dict.update({dump_dict_key:dump_dict_value})

	with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'wb') as dump_file:
		dill.dump(dump_dict, dump_file)

if __name__ == "__main__":
	prepare()
