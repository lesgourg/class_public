'''
This model has two input parameters:
    1) mass (between 50 and 100 GeV)
    2) the relative contribution of muons as primaries
'''

import numpy as np
import itertools
import dill
import os
import sys
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )

import DarkAges
from DarkAges.common import print_error, finalize, channel_dict
from DarkAges.model import annihilating_model, decaying_model
from DarkAges import redshift, options, transfer_functions, options

def run( *arguments ):
	##### In this block the external parameters in arguments[1:] (sys.argv[1:]) are read and translated into 
	## the parameters you like. In the given example case study the parameters are the DM mass at which the
	## spectra should be sampled and the relative contribbution of muons to the total spectrum.
	## If the input is not useful (mass outside the considered range, mixing not in [0,1], errors are raised
	#####
	if len(arguments) <3:
		print_error("There are too few arguments passed. I expected at least 2")
	sampling_mass = float(arguments[1])
	if sampling_mass < 50 or sampling_mass > 100:
		print_error("The mass-parameter sholud be in the range [50 GeV, 100 GeV]")
	mixing = float(arguments[2])
        #from warnings import warn
        #warn( "mixing = {:g}".format(mixing))
	if mixing < 0 or mixing > 1:
		print_error("The mixing-parameter sholud be in the range [0,1]")
	#####

	##### In this block we read in the the prepared 'interpolation'-objects
	#####
	model_dir = os.path.split(os.path.realpath(__file__))[0]
	model_name =  model_dir.split('/')[-1]

	with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'rb') as dump_file:
		dump_dict = dill.load(dump_file)
		spec_interp_muon = dump_dict.get('spec_interp_muon')
		spec_interp_bottom = dump_dict.get('spec_interp_bottom')
	####

	#### In this block the spectra are sampled and the total spectra are calculated
	#####
	temp_spec_muon = spec_interp_muon.__call__(sampling_mass)
	temp_spec_bottom = spec_interp_bottom.__call__(sampling_mass)

	total_spec = mixing * temp_spec_muon + (1.-mixing) * temp_spec_bottom
	#####

	##### In this block the 'model' object is created given the spectra (and history, decay time....)
	#####
	history = options.get('injection_history','annihilation')
	if history == 'decay':
		tdec = options.get('t_dec')
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
             f_functions[channel_dict['LowE']])
	#####

if __name__ == "__main__":
	run( *sys.argv )
