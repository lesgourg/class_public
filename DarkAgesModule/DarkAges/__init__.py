u"""
.. module:: DarkAges
   :synopsis: Module for the calculation of f(z) given a certain spectrum of particles
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

When loaded by :code:`import DarkAges` or by :code:`from DarkAges import ...`
"""

import numpy as np
np.seterr(all='ignore')
#Silence numpy when it comes to overflows etc.
import sys
import os

os.environ['DARKAGES_BASE'] = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.insert(0, os.environ['DARKAGES_BASE'])

logEnergies = None
redshift = None
transfer_functions = None
CosmoBackground = None

from .transfer import transfer, transfer_dump, transfer_load

DarkOptions = dict()

channel_dict = {
        'H-Ion': 0,
       'He-Ion': 1,
         'Ly-A': 2,
         'Heat': 3,
         'LowE': 4
}
"""Dictionary to translate between the order of deposition channels``
given by T. Slatyer and the order of channels used
in `CLASS <http://class-code.net>`_

+-------------------+-------+--------+-------+--------+------+
| index / column    | 0     |  1     | 2     | 3      | 4    |
+===================+=======+========+=======+========+======+
| **Slatyer order** | H-Ion | He-Ion | Ly-A  | Heat   | lowE |
+-------------------+-------+--------+-------+--------+------+
| **CLASS-order**   | Heat  | Ly-A   | H-Ion | He-Ion | lowE |
+-------------------+-------+--------+-------+--------+------+
"""

def print_info(message):
	u"""Covers messages for the informational use with #, such that the output
	is not caught by `CLASS <http://class-code.net>`_

	Parameters
	----------
	meassge : :obj:`str`
		Message to be printed in :obj:`stdout`
	"""

	print('#INFO: {0}'.format(message))

def print_warning(message):
	u"""Warning at critical behaviour of the code.

	Parameters
	----------
	meassge : :obj:`str`
		Message to be printed in :obj:`stderr` as :obj:`RuntimeWarning`
	"""

	from warnings import warn
	warn('\n\nWARNING: {0}\n'.format(message), RuntimeWarning )

class DarkAgesError(Exception):
	def __init__(self, message):
		Exception.__init__(self)
		self.message = '\n\n!!! ERROR !!!\n\n --> {} \n'.format(message)
		#print('\n\n ERROR:\n !!! {} !!!\n'.format(message))

	def __str__(self):
		return self.message


def get_background(key=None):
	if CosmoBackground is None:
		raise DarkAgesError('The background parameters are not set.')
	if key is None:
		return CosmoBackground
	else:
		if key in CosmoBackground.keys():
			return CosmoBackground.get(key)
		else:
			raise DarkAgesError('CosmoBackground has no key "{0}"'.format(key))

def set_background(H0 = 67.27, Om_M = 0.3156, Om_R = 8e-5):
	u"""Defines the background parameters :math:`H_0`, :math:`\\Omega_\\mathrm{matter}`
	and :math:`\\Omega_\\mathrm{radiation}`
	(Mostly for the use of :func:`H <DarkAges.common.H>`).

	Parameters
	----------
	H0 : :obj:`float`
		Todays Hubble parameter :math:`H_0`
		(*in units of* :math:`\\frac{\\mathrm{km}}{\\mathrm{Mpc\\,s}}`)
	Om_M : :obj:`float`
		Todays matter fraction :math:`\\Omega_\\mathrm{matter}`
	Om_R : :obj:`float`
		Todays radiation fraction :math:`\\Omega_\\mathrm{radiation}`

	Returns
	-------
	:obj:`tuple` of :obj:`float`
		tuple of :math:`H_0` (*in units of* :math:`\\frac{1}{\\mathrm{s}}`),
		:math:`\\Omega_\\mathrm{matter}` and :math:`\\Omega_\\mathrm{radiation}`
	"""

	_km_per_Mpc = 3.241e-20

	global CosmoBackground
	CosmoBackground = dict()
	CosmoBackground.update({'H0':H0*_km_per_Mpc,'Omega_m':Om_M,'Omega_r':Om_R})
	return

def get_redshift():
	if redshift is None:
		raise DarkAgesError('"redshift" was not initialized yet.')
	else:
		return redshift

def set_redshift(z):
	global redshift
	# As of now it makes only sense to use the redshift of the transfer functions
	# rather than a custom array, because the interpolation of the transfer functions
	# in redshift space should be done with caution and is not included in this
	# Version of the module. Therefore we only allow that this function can
	# be exectued once when the initial value of redshift is set.
	if redshift is None:
		redshift = z

def get_logEnergies():
	if logEnergies is None:
		raise DarkAgesError('"logEnergies" was not initialized yet.')
	else:
		#print "Someone asked for logEnergies (len = {:d})".format(logEnergies.__len__())
		return logEnergies

def set_logEnergies(logE):
	global logEnergies
	#if logEnergies is None:
	#	print "logEnergies set to its initial value (len = {:d})".format(logE.__len__())
	#else:
	#	print "logEnergies was updated (len = {:d})".format(logE.__len__())
	logEnergies = logE

def _transfer_init_and_dump(transfer_functions):
	for channel in channel_dict.keys():
		idx = channel_dict.get(channel)
		transfer_functions[idx] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch{:d}.dat'.format(idx+1)))
		transfer_dump(transfer_functions[idx], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch{:d}.obj'.format(idx+1)))
	return transfer_functions

def _transfer_load_from_dump(transfer_functions):
	for channel in channel_dict.keys():
		idx = channel_dict.get(channel)
		transfer_functions[idx] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch{:d}.obj'.format(idx+1)) )
	return transfer_functions

#################################

if transfer_functions is None:
	transfer_functions = np.empty(shape=5, dtype=transfer)

	transfer_is_initialized = True
	for i in xrange(5):
		transfer_is_initialized = transfer_is_initialized and os.path.isfile(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch{:d}.obj'.format(i+1)))

	if not transfer_is_initialized:
		print_info('The transfer seem not to be initialized. This will be done now. this may take a few seconds.')
		_transfer_init_and_dump(transfer_functions)
		print_info('The transfer functions are now initialized and loaded.\n')
	else:
		#print_info('The transfer functions are already initialized and loaded.\n')
		_transfer_load_from_dump(transfer_functions)
	del transfer_is_initialized
	del i

if logEnergies is None: set_logEnergies(transfer_functions[0].log10E[:])
if redshift is None: set_redshift(transfer_functions[0].z_deposited[:])
if CosmoBackground is None: set_background()

if 'DARKAGES_TOPLEVEL_PID' in os.environ:
	if os.getppid() == int(os.environ['DARKAGES_TOPLEVEL_PID']):
		import yaml
		top_level_random = int(os.environ['DARKAGES_TOPLEVEL_RANDN'])
		with open(os.path.join( os.environ['DARKAGES_BASE'], 'pid_{:d}_{:d}.yaml'.format(os.getppid(), top_level_random)), 'rb') as options_dumper:
			loaded_options = yaml.load(options_dumper)
			DarkOptions.update( loaded_options )
			set_background(DarkOptions.get('background_H0'),DarkOptions.get('background_Omega_m'),DarkOptions.get('background_Omega_r'))
			del loaded_options
			del top_level_random
