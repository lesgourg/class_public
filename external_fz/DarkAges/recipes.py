import numpy as np
from common import *
from model import model 

##### Functions related to executing a script-like file

def execute_script_file(ext_script_file, *arguments):
	import subprocess
	command = ['python']
	command.append(ext_script_file)
	if arguments:
		for arg in arguments[0]:
			command.append(arg)
	print_info('running script-file: "{}"'.format(ext_script_file))
	retcode = subprocess.call(command)
	if retcode != 0:
		print_error('Failed to execute the script-file: "{}"'.format(ext_script_file))

##### Functions related to loading a model from a file contiaining the input spectra (and mass)

def loading_from_specfile(fname, transfer_functions, logEnergies, redshift, history):
	model_from_file = load_from_spectrum(fname, logEnergies, hist=history)
	try:
		assert len(channel_dict) == len(transfer_functions)
	except AssertionError:
		print_error('The number of "transfer" instances ({.d}) and the number of channels ({.d}) do not match'.format(len(transfer_functions),len(channel_dict)))
	f_function = np.zeros( shape=(len(channel_dict),len(redshift)), dtype=np.float64 )
	for channel in channel_dict:
		idx = channel_dict[channel]
		f_function[idx,:] = model_from_file.calc_f(transfer_functions[idx])[-1]

	finalize(redshift, 
       		 f_function[channel_dict['Heat']], 
       		 f_function[channel_dict['Ly-A']], 
       		 f_function[channel_dict['H-Ion']],
       		 f_function[channel_dict['He-Ion']],        
			 f_function[channel_dict['LowE']])

def load_from_spectrum(fname,  logEnergies, hist='annihilation'):
	spec_data = np.genfromtxt(fname, unpack=True, usecols=(0,1,2,3,4), skip_header=1, dtype=np.float64)
	masses = np.unique( spec_data[0] )
	if masses.shape[-1] > 1:
		example_model = np.empty(shape = masse.shape, dtype=model)
		for idx in xrange(mass.shape[-1]):
			mass = masses[idx]
			mass_mask = (masses == mass) 
			spec_el, spec_ph, spec_oth = sample_spectrum(spec_data[2,mass_mask], spec_data[3,mass_mask], spec_data[4,mass_mask], spec_data[1,mass_mask], mass, logEnergies)
			#example_model[idx] = model(spec_el, spec_ph, spec_oth, 1e9*mass)
			example_model[idx] = model(spec_el, spec_ph, spec_oth, 1e9*mass, t_dec = 1e18, history=hist)
	else:
		mass = masses[0]
		spec_el, spec_ph, spec_oth = sample_spectrum(spec_data[2], spec_data[3], spec_data[4], spec_data[1], mass, logEnergies)
		#example_model = model(spec_el, spec_ph, spec_oth, 1e9*mass)
		example_model = model(spec_el, spec_ph, spec_oth, 1e9*mass, t_dec = 1e18, history=hist)

	return example_model

##### Functions related to running a preprocessed model (or defining it, if it does not exist)

def access_model(model_name, *arguments):
	model_dir = os.path.join(os.environ['DARKAGES_BASE'], 'models/{}'.format(model_name))
	if os.path.isfile( os.path.join(model_dir, '{}.obj'.format(model_name)) ):
		if arguments:	
			run_model(model_dir, arguments[0])
		else:
			run_model(model_dir)
	else:
		prepare_model(model_dir)

def prepare_model(model_dir):
	import subprocess
	command = ['python']
	file_to_run = os.path.join(model_dir,'prepare.py')
	command.append(file_to_run)
	print_info('Preparing_the model: running script: "{}"'.format(file_to_run))
	retcode = subprocess.call(command)
	if retcode != 0:
		print_error('Failed to prepare the model. Error in the execution of "{}"'.format(file_to_run))
	else:
		print_info('Finished preparing the model. It is now ready to use. Please rerun your command.')

def run_model(model_dir, *arguments):
	import subprocess
	command = ['python']
	file_to_run = os.path.join(model_dir,'run.py')
	command.append(file_to_run)
	if arguments:
		for arg in arguments[0]:
			command.append(arg)
	print_info('running script-file: "{}"'.format(file_to_run))
	retcode = subprocess.call(command)
	if retcode != 0:
		print_error('Failed to execute the script-file: "{}"'.format(file_to_run))

	
