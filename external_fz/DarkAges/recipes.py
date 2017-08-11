import numpy as np
import os
import sys
from .common import channel_dict, finalize, sample_spectrum, print_info, print_error, print_warning
from .__init__ import redshift, logEnergies, transfer_functions, options
from .model import model, annihilating_model, decaying_model, evaporating_model
from .interpolator import logInterpolator, NDlogInterpolator

##### Functions related to executing a script-like file

def execute_script_file(ext_script_file, *arguments):
	import subprocess
	command = ['{0}'.format(sys.executable)]
	#command.append('-OO')
	command.append(ext_script_file)
	if arguments:
		for arg in arguments[0]:
			command.append(arg)
	print_info('running script-file: "{0}"'.format(ext_script_file))
	retcode = subprocess.call(command)
	if retcode != 0:
		print_error('Failed to execute the script-file: "{0}"'.format(ext_script_file))

##### Functions related to loading a model from a file containing the input spectra (and mass)

def evaporating_PBH( PBH_mass_ini, transfer_functions, logEnergies, redshift , merge_ion = False):
	model_from_file = evaporating_model(PBH_mass_ini)
	f_function = np.zeros( shape=(len(channel_dict),len(redshift)), dtype=np.float64 )
	for channel in channel_dict:
		idx = channel_dict[channel]
		f_function[idx,:] = model_from_file.calc_f(transfer_functions[idx])[-1]

	if merge_ion:
		f_ion_lya = f_function[channel_dict['Ly-A'],:]
		f_ion_H = f_function[channel_dict['H-Ion'],:]
		f_ion_He = f_function[channel_dict['He-Ion'],:]
		f_function[channel_dict['H-Ion'],:] = np.sum(np.asarray([f_ion_lya, f_ion_H, f_ion_He]), axis=0)
		f_function[channel_dict['Ly-A'],:] = 0.
		f_function[channel_dict['He-Ion'],:] = 0.

	finalize(redshift,
			 f_function[channel_dict['Heat']],
			 f_function[channel_dict['Ly-A']],
			 f_function[channel_dict['H-Ion']],
			 f_function[channel_dict['He-Ion']],
			 f_function[channel_dict['LowE']],
			 **options)


def loading_from_specfiles(fnames, transfer_functions, logEnergies, redshift, mass, t_dec, hist='annihilation', branchings=[1.], **options):
	branchings = np.asarray(branchings)
	spectra = np.empty(shape=(3,len(logEnergies),len(fnames)), dtype=np.float64)
	for idx, fname in enumerate(fnames):
		spec_interpolator = load_from_spectrum(fname, logEnergies, **options)
		lower = spec_interpolator.get_lower()
		upper = spec_interpolator.get_upper()
		if mass < lower or mass > upper:
			print_warning('The spectra-file >>{:s}<< contains only spectra in the mass range [{:.2g}, {:.2g}]. Hence the spectrum you asked for (mass: {:.2g}) cannot be deduced. Return zeros'.format(fname, lower, upper, mass))
			spectra[:,:,idx] = np.zeros(shape=(3,len(logEnergies)), dtype=np.float64)
		else:
			spectra[:,:,idx] = spec_interpolator.__call__(mass)
	try:
		assert spectra.shape[-1] == branchings.shape[-1]
	except AssertionError:
		print_error('The number of spectra ({:d}) and the number of provided branching ratios ({:d}) do not match'.format(spectra.shape[-1],branchings.shape[-1]))
	tot_spec = np.tensordot(spectra, branchings, axes=(2,0))
	if hist == 'decay':
		#model_from_file = model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass, t_dec = t_dec, history=hist)
		model_from_file = decaying_model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass, t_dec)
	elif hist == 'annihilation':
		#model_from_file = model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass, history=hist)
		model_from_file = annihilating_model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass)
	else:
		print_error('The method >> {:s} << cannot deal with the injection history >> {:s} <<'.format(loading_from_specfiles.func_name, hist))
	try:
		assert len(channel_dict) == len(transfer_functions)
	except AssertionError:
		print_error('The number of "transfer" instances ({:d}) and the number of channels ({:d}) do not match'.format(len(transfer_functions),len(channel_dict)))
	f_function = np.zeros( shape=(len(channel_dict),len(redshift)), dtype=np.float64 )
	for channel in channel_dict:
		idx = channel_dict[channel]
		f_function[idx,:] = model_from_file.calc_f(transfer_functions[idx])[-1]

	finalize(redshift,
       		 f_function[channel_dict['Heat']],
       		 f_function[channel_dict['Ly-A']],
       		 f_function[channel_dict['H-Ion']],
       		 f_function[channel_dict['He-Ion']],
			 f_function[channel_dict['LowE']],
             **options)

def load_from_spectrum(fname, logEnergies, **options):
	cols_to_use = options.get('spectra_cols',(0,1,2,3,4))
	try:
		assert len(cols_to_use) == 5
	except AssertionError:
		print_error('There is something wrong with the number of columns in the spectra-file. I expected 5 columns but I got only {:d}'.format(len(cols_to_use)))
	spec_data = np.genfromtxt(fname, unpack=True, usecols=(0,1,2,3,4), skip_header=1, dtype=np.float64)
	masses = np.unique(spec_data[0,:])
	temp_spec = np.empty(shape=(len(masses),3,len(logEnergies)), dtype=np.float64)
	for idx, mass in enumerate(masses):
		mass_mask = np.absolute(spec_data[0] - mass) <= 1e-5
		temp_spec[idx,:,:] = sample_spectrum(spec_data[2,mass_mask], spec_data[3,mass_mask], spec_data[4,mass_mask], spec_data[1,mass_mask], mass, logEnergies, **options)

	return NDlogInterpolator(masses, temp_spec, exponent=1)

def load_from_spectrum2(fnames, logEnergies, mass, t_dec, hist='annihilation', branchings=np.asarray([1.])):
	if type(fnames) is not list and type(fnames) is not np.ndarray:
		fnames = np.asarray([fnames])
	else:
		fnames = np.asarray(fnames)
	temp_spec = np.empty(shape=(3,len(fnames),len(logEnergies)))
	for idx, fname in enumerate(fnames):
		spec_data = np.genfromtxt(fname, unpack=True, usecols=(0,5,6,7,8), skip_header=1, dtype=np.float64)
		mass_mask = np.absolute(spec_data[0] - mass) <= 1e-5
		temp_spec[:,idx,:] = sample_spectrum(spec_data[2,mass_mask], spec_data[3,mass_mask], spec_data[4,mass_mask], spec_data[1,mass_mask], mass, logEnergies)
	spec_el, spec_ph, spec_oth = np.tensordot(temp_spec, branchings, axes=(1,0))
	if hist == 'decay':
		example_model = model(spec_el, spec_ph, spec_oth, 1e9*mass, t_dec = t_dec, history='decay')
	else:
		example_model = model(spec_el, spec_ph, spec_oth, 1e9*mass, history=hist)

	return example_model

##### Functions related to running a preprocessed model (or defining it, if it does not exist)

def access_model(model_name, force_rebuild = False, *arguments):
	model_dir = os.path.join(os.environ['DARKAGES_BASE'], 'models/{0}'.format(model_name))
	sys.path.insert(0,model_dir)
	if os.path.isfile( os.path.join(model_dir, '{0}.obj'.format(model_name)) ) and not force_rebuild:
		run_model(model_dir, *arguments)
	else:
		prepare_model(model_dir)

def prepare_model(model_dir):
	file_to_run = os.path.join(model_dir,'prepare.py')
	_temp = __import__('prepare', globals(), locals(), [], -1)
	_prepare = _temp.prepare
	print_info('Preparing_the model: running script: "{0}"'.format(file_to_run))
	_prepare()
	print_info('Finished preparing the model. It is now ready to use. Please rerun your command.')

def run_model(model_dir, *arguments):
	_cmnd = []
	file_to_run = os.path.join(model_dir,'run.py')
	_cmnd.append(file_to_run)
	if arguments:
		for arg in arguments[0]:
			_cmnd.append(arg)
	_temp = __import__('run', globals(), locals(), [], -1)
	_run = _temp.run
	print_info('running script-file: "{0}"'.format(file_to_run))
	_run(*_cmnd)
	
