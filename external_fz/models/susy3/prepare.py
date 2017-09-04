import itertools
import dill
import numpy as np
import sys
import os
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges import logEnergies, DarkOptions
from DarkAges.recipes import load_from_spectrum

#####
def prepare():
	model_dir = os.path.split(os.path.realpath(__file__))[0]
	model_name =  model_dir.split('/')[-1]

	primaries = np.array(['onlyZ','WandZ'], dtype=np.dtype('a16'))

	dump_dict = dict()

	for idx_prim, primary in enumerate(primaries):
		fname = os.path.join(model_dir, 'data/susy3_right_{:s}.dat'.format(primary))
		dump_dict_key = 'spec_interp_right_{:s}'.format(primary)
		dump_dict_value = load_from_spectrum(fname, logEnergies, **DarkOptions)
		dump_dict.update({dump_dict_key:dump_dict_value})

	with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'wb') as dump_file:
		dill.dump(dump_dict, dump_file)
#####
if __name__ == "__main__":
	prepare()
