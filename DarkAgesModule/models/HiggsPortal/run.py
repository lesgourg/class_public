import numpy as np
import dill
import os
import sys
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges.common import finalize
from DarkAges.model import annihilating_model, decaying_model
from DarkAges import redshift, DarkAgesError, transfer_functions, channel_dict

model_dir = os.path.split(os.path.realpath(__file__))[0]
model_name =  model_dir.split('/')[-1]

def run( *arguments, **DarkOptions ):
	if len(arguments) <2:
		raise DarkAgesError("There are too few arguments passed. I expected at least 1")
	sampling_mass = float(arguments[1])
	if sampling_mass < 5 or sampling_mass > 5e3:
		raise DarkAgesError("The mass-parameter sholud be in the range [5 GeV, 5 TeV]")

	with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'rb') as dump_file:
		dump_dict = dill.load(dump_file)

	primaries = dump_dict.get('channels')
	tmp_shape = dump_dict.get('spec_interp_{:s}'.format(primaries[0]))._shape
	total_spec = np.zeros(shape=tmp_shape, dtype=np.float64)
	tot_frac = 0.
	for primary in primaries:
		fraction = dump_dict.get('frac_interp_{:s}'.format(primary)).__call__(sampling_mass)
		#fraction = dump_dict.get('frac_interp_{:s}_loglin'.format(primary)).__call__(sampling_mass)
		spectra =  dump_dict.get('spec_interp_{:s}'.format(primary)).__call__(sampling_mass)
		if np.any(spectra > 0.0):
			total_spec += fraction * spectra
			tot_frac += fraction
	if tot_frac > 0.:
		total_spec /= tot_frac
	else:
		total_spec = np.zeros_like(total_spec)

	history = DarkOptions.get('injection_history','annihilation')
	if history == 'decay':
		tdec = DarkOptions.get('t_dec')
		full_model = decaying_model(total_spec[0], total_spec[1], total_spec[2], 1e9*sampling_mass, tdec)
	else:
		full_model = annihilating_model(total_spec[0], total_spec[1], total_spec[2], 1e9*sampling_mass)
	#####

	##### To finish the calculation, calculate f(z) for each deposition channel
	## and print the CLASS-output with the finalize()-method
	#####
	f_functions = np.zeros((len(channel_dict),len(redshift)))
	for channel in channel_dict:
		idx = channel_dict[channel]
		f_functions[idx,:] = full_model.calc_f(transfer_functions[idx])[-1]

	finalize(redshift,
             f_functions[channel_dict['Heat']],
             f_functions[channel_dict['Ly-A']],
             f_functions[channel_dict['H-Ion']],
             f_functions[channel_dict['He-Ion']],
             f_functions[channel_dict['LowE']],
			 **DarkOptions)
	#####

if __name__ == "__main__":
	from DarkAges import DarkOptions as milk
	run( *sys.argv, **milk)
