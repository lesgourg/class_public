import itertools
import dill
import numpy as np
import sys
import os
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges import get_logEnergies
from DarkAges.recipes import load_from_spectrum

def prepare():
	##### In this Block, the model is introduced
	## e.g. Introducing the folder where the data is
	## or introducing the names of the primaries to consider,
	## so that there spectra can be found and read
	#####
	model_dir = os.path.split(os.path.realpath(__file__))[0]
	model_name =  model_dir.split('/')[-1]

	primaries = np.array(['muon','bottom'], dtype=np.dtype('U32'))
	#####

	##### The idea of the 'model'-mode is to do the interpolation of the spectra once
	## and then save the 'interpolator', that the spectrum can then easily be sampled
	## by loading the object and run the .__call__ method.
	## To this end, we store the interpolators in a dictionary and save this
	## dictionary with the dill-module.
	## To get the energies at which we need to sample the spectra we make
	## use of the get_logEnergies-method which imports the global logEnergies array
	#####
	dump_dict = dict()
	logEnergies = get_logEnergies()

	for idx_prim, primary in enumerate(primaries):
		fname = os.path.join(model_dir, 'data/{:s}_50-100_spectrum.dat'.format(primary))
		dump_dict_key = 'spec_interp_{:s}'.format(primary)
		dump_dict_value = load_from_spectrum(fname, logEnergies)
		dump_dict.update({dump_dict_key:dump_dict_value})

	with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'wb') as dump_file:
		dill.dump(dump_dict, dump_file)
	#####

if __name__ == "__main__":
	prepare()
