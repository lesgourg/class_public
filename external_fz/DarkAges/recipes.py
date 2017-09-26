u"""
.. module:: recipes
   :synopsis: Definition of "Standard recipes" to calculate f(z)
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

Contains the definitions of wrapper to the most used scenarios to calculate
:math:`f(z)`. This methods are mainly imported by the caller (:code:`/DARK_AGES_ROOT/bin/DarkAges`)
and automatically initiated depending on the parsed inputs.

By now there are three possible scenarios and the following wrapping methods:

- Executing an external Python-script which uses methods of the DarkAges-package:

	* :meth:`execute_script_file`

- Reading spectra from external tables or calculating an analytical spectrum
  and setting up an instance of :class:`model <DarkAges.model.model>` and
  calculate :math:`f(z)`:

    * :meth:`accreting_PBH`: Accretion of a primordial black hole
	  The spectrum is calculated according to the mass-loss
	  (see :class:`evaporator <DarkAges.evaporator>`), the
	  :class:`model <DarkAges.model.model>` with the appropriate injection
	  history is initialized and :math:`f(z)` is calculated.
	* :meth:`evaporating_PBH`: Evaporation of a primordial black hole
	  The spectrum is calculated according to the mass-loss
	  (see :class:`evaporator <DarkAges.evaporator>`), the
	  :class:`model <DarkAges.model.model>` with the appropriate injection
	  history is initialized and :math:`f(z)` is calculated.
	* :meth:`loading_from_specfiles`: Spectra are read from external tables,
	  and an instance of :class:`model <DarkAges.model.model>` with the
	  appropriate injection history is initialized and :math:`f(z)` is calculated.
	* :meth:`load_from_spectrum` :  Spectra are read from external tables,
	  an instance of :class:`NDlogInterpolator <DarkAges.interpolator.NDlogInterpolator>`
	  is initialized and returned

- Calling (or initializing if called the first time) a standardized Python-script
  in the "models"-folder for a custom model (See LINK TO THE SECTION)

	* :meth:`access_model`: Access the model in a given subfolder of "models".
	  Either call the :func:`run` or :func:`prepare`-routine.

"""

import numpy as np
import os
import sys
from .common import channel_dict, finalize, sample_spectrum, print_info, print_warning
#from .__init__ import DarkOptions as options
from .__init__ import redshift, logEnergies, transfer_functions, DarkAgesError
from .model import model, annihilating_model, decaying_model, evaporating_model, annihilating_halos_model, accreting_model
from .interpolator import logInterpolator, NDlogInterpolator

##### Functions related to executing a script-like file

def execute_script_file(ext_script_file, *arguments):
	u"""This method forks a childprocess which runs a new python interpreter
	on a Python-script which uses the DarkAges-package (or at least parts of it).

	To do so the :meth:`call`-method of the :class:`subprocess`-package is used.

	The main purpose of this method is to run the script with updated background
	parameters given by :code:`/DARK_AGES_ROOT/bin/DarkAges --cosmo_bg .....`

	.. warning:: We discourage discourage you to heavily run the code in this
	   configuration like for parameter extraction with MontePython, since it
	   is very IO-intensive with updating of the background parameter and further
	   options via a temporary file.
	   This leads to an bottleneck in the execution. We advise you to consider
	   the possibility to rewrite your routine in the script in a way that it can
	   be understood as a part of the "models"-routines (See :ref:`using_the_own_model`).

	Parameters
	----------
	ext_script_file : :obj:`str`
		Filename (absolute or relative) to the Python script to execute.
	*arguments : :obj:`tuple`
		Additional parameters to the script (Are usually passed as
		:code:`sys.argv[1:]`)
	"""

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
		raise DarkAgesError('Failed to execute the script-file: "{0}"'.format(ext_script_file))

##### Functions related to loading a model from a file containing the input spectra (and mass)
def accreting_PBH( PBH_mass, recipe, transfer_functions, logEnergies, redshift , merge_ion = False, **DarkOptions):
    u"""Wrapper for the calculation of :math:`f_c (z)` for a evaporating primordial
    black hole (PBH) with a given initial mass :code:`PBH_mass_ini` and prionts
    the table of them for all five deposition channels

    Optionally the channels the channels 'Ly-A excitation', 'helium ionization',
    and 'hydrogen ionization' are merged into the hydrogen ionization channel and
    the other two channels are left blank.

    Parameters
    ----------
    PBH_mass : :obj:`float`
    	Mass of the primordial black hole (*in units of* :math:`M_\odot`)
    recipe : :obj:`string`
        Recipe setting the luminosity and the rate of the accretion (`spherical_accretion` taken from 1612.05644 and `disk_accretion` from 1707.04206)
    logEnergies : :obj:`array-like`, optional
    	Array (:code:`shape = (l)`) of the logarithms of the kinetic energies of the particles
    	(*in units of* :math:`\\mathrm{eV}`) to the base 10.
    	If not specified, the standard array provided by
    	:class:`the initializer <DarkAges.__init__>` is taken.
    redshift : :obj:`array-like`, optional
    	Array (:code:`shape = (k)`) with the values of :math:`z+1`. Used for
    	the calculation of the double-differential spectra.
    	If not specified, the standard array provided by
    	:class:`the initializer <DarkAges.__init__>` is taken.

    """

    model_from_file = accreting_model(PBH_mass,recipe)
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
             **DarkOptions)


def evaporating_PBH( PBH_mass_ini, transfer_functions, logEnergies, redshift , merge_ion = False, **DarkOptions):
	u"""Wrapper for the calculation of :math:`f_c (z)` for a evaporating primordial
	black hole (PBH) with a given initial mass :code:`PBH_mass_ini` and prionts
	the table of them for all five deposition channels

	Optionally the channels the channels 'Ly-A excitation', 'helium ionization',
	and 'hydrogen ionization' are merged into the hydrogen ionization channel and
	the other two channels are left blank.

	Parameters
	----------
	PBH_mass_ini : :obj:`float`
		Initial mass of the primordial black hole (*in units of* :math:`\\mathrm{g}`)
	transfer_functions : :obj:`class`
		Array of initialized instance of :class:`transfer <DarkAges.transfer.transfer>`
		for the five deposition channels in question in the order imposed by
		T. Slatyer (c.f. :meth:`channel_dict <DarkAges.common.channel_dict>`)
	logEnergies : :obj:`array-like`
		Array (:code:`shape = (l)`) with the values of the logarithm to the base 10
		of the kinetic energy of the particles
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) with the values of :math:`z+1` of the redshift
		of injection
	merge_ion : :obj:`bool`, *optional*
		Flag to specify if the channels refering to ionization should be taken
		as one channel and the other channels should be left blank.
		If not given the default of :code:`merge_ion = False` is taken.
	"""

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
			 **DarkOptions)


def loading_from_specfiles(fnames, transfer_functions, logEnergies, redshift, mass, t_dec, zh=1.,fh=0., hist='annihilation', branchings=[1.], **DarkOptions):
	u"""Wrapper to calculate :math:`f(z)` and print the table for all five deposition channels
	from spectra tabulated in files for a given injection history.

	If more than one spectrum is given, an information about their relative weighting needs to be
	given by :code:`branching`.

	Parameters
	----------
	fnames : :obj:`str`
		Filename (or array of filenames) which include the table of the spectra.
		Can either be relative or absolute.
	transfer_functions : :obj:`class`
		Array of initialized instance of :class:`transfer <DarkAges.transfer.transfer>`
		for the five deposition channels in question in the order imposed by
		T. Slatyer (c.f. :meth:`channel_dict <DarkAges.common.channel_dict>`)
	logEnergies : :obj:`array-like`
		Array (:code:`shape = (l)`) with the values of the logarithm to the base 10
		of the kinetic energy of the particles
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) with the values of :math:`z+1` of the redshift
		of injection
	mass : :obj:`float`
		Mass of the DM candidate (*in units of* :math:`\\mathrm{GeV}`)
	t_dec : :obj:`float`
		Lifetime of the DM candidate for a decaying species (*in units of* :math:`\\mathrm{s}`)
		Mandatory for :code:`hist='decay'` (Will be ignored for
		:code:`hist='annihilation'`. In taht case you can set it to infinity)
	hist : :obj:`str`, *optional*
		String with the energy injection history to consider. Valid options are:
		:code:`'annihilation'`, and  :code:`'decay'`. If not given the default
		value :code:`hist='annihilation'` will be taken
	branchings : :obj:`array-like`
		Array of the relative contributions of each of the spectra specified in
		:code:`fnames`. Need to add upp to 1. and need to have the same number
		of entries as :code:`fnames`.
		If :code:`fnames` has on one entry, the default value :code:`branchings=[1.]`
		is taken

	Raises
	------
	DarkAgesError
		if the entries in :code:`branchings` do not add up to one and/or the
		number of entries in :code:`branchings` is not consistent with the number of
		entries in :code:`fnames`
	"""

	branchings = np.asarray(branchings)
	spectra = np.empty(shape=(3,len(logEnergies),len(fnames)), dtype=np.float64)
	if fnames != ['Dirac'] and fnames != ['dirac']:
		for idx, fname in enumerate(fnames):
			spec_interpolator = load_from_spectrum(fname, logEnergies, **DarkOptions)
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
           		raise DarkAgesError('The number of spectra ({:d}) and the number of provided branching ratios ({:d}) do not match'.format(spectra.shape[-1],branchings.shape[-1]))
        	tot_spec = np.tensordot(spectra, branchings, axes=(2,0))
	elif fnames == ['Dirac'] or fnames == ['dirac']:
		print 'here'

	if hist == 'decay':
		model_from_file = decaying_model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass, t_dec)
	elif hist == 'annihilation':
		model_from_file = annihilating_model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass)
	elif hist == 'annihilation_halos':
		model_from_file = annihilating_halos_model(tot_spec[0], tot_spec[1], tot_spec[2], 1e9*mass,zh,fh)
	else:
		raise DarkAgesError('The method >> {:s} << cannot deal with the injection history >> {:s} <<'.format(loading_from_specfiles.func_name, hist))
	try:
		assert len(channel_dict) == len(transfer_functions)
	except AssertionError:
		raise DarkAgesError('The number of "transfer" instances ({:d}) and the number of channels ({:d}) do not match'.format(len(transfer_functions),len(channel_dict)))
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
             **DarkOptions)

def load_from_spectrum(fname, logEnergies, **DarkOptions):
	u"""Wrapper to return the interpolated spectra of electrons and positrons,
	photons, and other particles as an initialized instance of the
	:class:`NDlogInterpolator <DarkAges.interpolator.NDlogInterpolator>`-class
	(with :code:`shape = (3,l)` for every of the three types of particles and every
	enrgy given in :code:`logEnergies`) given the filename of the table of
	spectra to read.

	The tables to read are expected to have the following shape

		+------+--------+--------------+--------------+---------------+
		| mass | log10E | dN/dE (elec) | dN/dE (phot) | dN/dE (other) |
		+------+--------+--------------+--------------+---------------+

	where the Energy is per default expected to be in units of :math:`\\mathrm{GeV}`.
	If your table differs in terms of that the energy scale is different from
	:math:`\\mathrm{GeV}`, that the spectra are given in terms of
	:math:`E \\cdot \\frac{\\mathrm{d}N}{\\mathrm{d}E}` rather than
	:math:`\\frac{\\mathrm{d}N}{\\mathrm{d}E}` or that the energies are given by
	their values and not their logarithm, please use the respective keyword arguments
	in :code:`**options` (see :meth:`sample_spectrum <DarkAges.common.sample_spectrum>`)

	It is expected that the file contains the tables for more than one
	DM-mass.

	For every single value of the DM-mass (the first column) of the table
	the spectra (columns 3-5) are interploated and sampled at the energies
	specified by the input :code:`logEnergies`. Afterwards the interpolations
	over the mass-array (first column) are initialized for every energy
	in the :code:`logEnergies`-array.

	Parameters
	----------
	fname : :obj:`str`
		Filename (relative or absolute) of the tabulated spectra.
	logEnergies : :obj:`array-like`
		Array (:code:`shape = (l)`) of the logarithms of the particles
		kinetic energies (*in units of* :math:`\\mathrm{eV}`) to the base 10
	**options : :obj:`Keyword arguments`
		Further options

	Returns
	-------
	:obj:`class`
		Initialized instance of the :class:`NDlogInterpolator <DarkAges.interpolator.NDlogInterpolator>`-class
		(:code:`shape = (3,l)`)

	Raises
	------
	DarkAgesError
		if the table has not exactly five columns.
	"""

	cols_to_use = DarkOptions.get('spectra_cols',(0,1,2,3,4))
	try:
		assert len(cols_to_use) == 5
	except AssertionError:
		raise DarkAgesError('There is something wrong with the number of columns in the spectra-file. I expected 5 columns but I got only {:d}'.format(len(cols_to_use)))
	spec_data = np.genfromtxt(fname, unpack=True, usecols=cols_to_use, skip_header=1, dtype=np.float64)
	masses = np.unique(spec_data[0,:])
	temp_spec = np.empty(shape=(len(masses),3,len(logEnergies)), dtype=np.float64)
	for idx, mass in enumerate(masses):
		mass_mask = np.absolute(spec_data[0] - mass) <= 1e-5
		temp_spec[idx,:,:] = sample_spectrum(spec_data[2,mass_mask], spec_data[3,mass_mask], spec_data[4,mass_mask], spec_data[1,mass_mask], mass, logEnergies, **DarkOptions)

	return NDlogInterpolator(masses, temp_spec, exponent=1)

##### Functions related to running a preprocessed model (or defining it, if it does not exist)

def access_model(model_name, force_rebuild = False, *arguments, **DarkOptions):
	u"""Wrapper to access the :meth:`run` and :meth:`prepare`-methods of
	the standardized python scripts in the :code:`models/...`-subfolders.

	For further information please read the section :ref:`using_the_own_model`

	Parameters
	----------
	model_name : :obj:`str`
		Name of the model (i.e. name of the :code:`models/...`-subfolder to search
		in)
	force_rebuild : :obj:`bool`, *optional*
		Flag to specify wether the :meth:`prepare`-method of the model should be
		rerun. (As long as 'model_name'.obj exists in the :code:`models/'model_name'`
		-subfolder the code assumes the model to be prepared and executes the
		:meth:`run`-method). The default value is set to :code:`False`
	*arguments : :obj:`tuple`
		Additional parameters to the script (Are usually passed as
		:code:`sys.argv[1:]`)
	"""

	model_dir = os.path.join(os.environ['DARKAGES_BASE'], 'models/{0}'.format(model_name))
	if not os.path.isdir(model_dir):
		raise DarkAgesError('"{0}" seems not to be a proper model. Could not find the respective directory.'.format(model_name))
	sys.path.insert(0,model_dir)
	if os.path.isfile( os.path.join(model_dir, '{0}.obj'.format(model_name)) ) and not force_rebuild:
		_run_model(model_dir, *arguments, **DarkOptions)
	else:
		_prepare_model(model_dir)

def _prepare_model(model_dir):
	_found = False
	for fname in os.listdir(model_dir):
		if fname.find('.py') != -1:
			try:
				if fname.find('pyc') != -1:
					_module = fname[:-4]
				else:
					_module = fname[:-3]
				_temp = __import__(_module, globals(), locals(), [], -1)
				_prepare = _temp.__getattribute__('prepare')
				file_to_run = os.path.join(model_dir,os.path.join(model_dir,fname))
				_found = True # If module with the "prepare"-method could be loaded exit the for loop
				break
			except AttributeError: # Proceed if the imported module has no function "prepare"
				pass
			except ImportError: # Proceed if the module cannot be imported
				pass
	if not _found: # If the loop finished and nothing was found raise an Error
		raise DarkAgesError('Could not find the method "prepare" in the directory of the model')
	print_info('Preparing_the model: executing prepare()-method found in "{0}".'.format(file_to_run))
	_prepare()
	print_info('Finished preparing the model. It is now ready to use. Please rerun your command.')

def _run_model(model_dir, *arguments, **DarkOptions):
	_cmnd = []
	_found = False
	for fname in os.listdir(model_dir):
		if fname.find('.py') != -1:
			try:
				if fname.find('pyc') != -1:
					_module = fname[:-4]
				else:
					_module = fname[:-3]
				_temp = __import__(_module, globals(), locals(), [], -1)
				_run = _temp.__getattribute__('run')
				file_to_run = os.path.join(model_dir,os.path.join(model_dir,_module,'.py'))
				_found = True
				break
			except AttributeError:
				pass
			except ImportError:
				pass
	if not _found:
		raise DarkAgesError('Could not find the method "run" in the directory of the model')
	_cmnd.append(file_to_run)
	if arguments:
		for arg in arguments[0]:
			_cmnd.append(arg)
	print_info('Executing run()-method found in "{0}".'.format(file_to_run))
	_run(*_cmnd, **DarkOptions)
