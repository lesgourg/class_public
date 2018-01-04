'''
This model has two input parameters:
    1) mass (between 50 and 500 GeV)
    2) the finalstate composition. Here: The allowed third particle (0 for 'onlyZ' and 1 for 'WandZ')
'''

import numpy as np
import itertools
import dill
import os
import sys
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges.common import finalize, channel_dict
from DarkAges.model import annihilating_model, decaying_model
from DarkAges import redshift, transfer_functions, DarkAgesError

#####
def run( *arguments, **DarkOptions ):
	if len(arguments) <3:
		raise DarkAgesError("There are too few arguments passed. I expected at least 2")
	#sampling_mass = float(arguments[1])
	sampling_mass = 5*(10**float(arguments[1]))
	if sampling_mass < 50 or sampling_mass > 500:
		raise DarkAgesError("The mass-parameter shsould be in the range [50 GeV, 500 GeV]")
	idx_primary = int(arguments[2])
	if idx_primary != 1 and idx_primary != 0:
		raise DarkAgesError("The mode '{:d}' is not recognized. Please enter a valid mode. The choices are: 0 for 'onlyZ', and 1 for 'WandZ'.".format(idx_primary))

	model_dir = os.path.split(os.path.realpath(__file__))[0]
	model_name =  model_dir.split('/')[-1]

	with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'rb') as dump_file:
		dump_dict = dill.load(dump_file)
		if idx_primary == 0:
			spec_interp = dump_dict.get('spec_interp_right_onlyZ')
		elif idx_primary == 1:
			spec_interp = dump_dict.get('spec_interp_right_WandZ')

	total_spec = spec_interp.__call__(sampling_mass)
	#print total_spec

	history = DarkOptions.get('injection_history','annihilation')
	if history == 'decay':
		tdec = DarkOptions.get('t_dec')
		full_model = decaying_model(total_spec[0], total_spec[1], total_spec[2], 1e9*sampling_mass, tdec)
	else:
		full_model = annihilating_model(total_spec[0], total_spec[1], total_spec[2], 1e9*sampling_mass)

	f_functions = np.zeros((len(channel_dict),len(redshift)))
	for channel in channel_dict:
		idx = channel_dict[channel]
		f_functions[idx,:] = full_model.calc_f(transfer_functions[idx])[-1]

	finalize(redshift,
             f_functions[channel_dict['Heat']],
             f_functions[channel_dict['Ly-A']],
             f_functions[channel_dict['H-Ion']],
             f_functions[channel_dict['He-Ion']],
             f_functions[channel_dict['LowE']])
#####
if __name__ == "__main__":
	from DarkAges import DarkOptions as milk
	run( *sys.argv, **milk )
