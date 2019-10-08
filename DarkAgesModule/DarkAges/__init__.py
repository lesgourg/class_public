u"""
.. module:: DarkAges
   :synopsis: Module for the calculation of f(z) given a certain spectrum of particles
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

When loaded by :code:`import DarkAges` or by :code:`from DarkAges import ...`
the default values of global variables used by the methods of the DarkAges-module
are set.

In particular, these variables are

*  :code:`CosmoBackground`: Dictionary containing the cosmological parameters
   :math:`H_0`, :math:`\\Omega_\\mathrm{m}` and :math:`\\Omega_\\mathrm{r}`,
   which can be acessed with the setter and getter functions
   :meth:`set_background` and :meth:`set_background`.
*  :code:`transfer_functions`: Array of initialized instances of the
   :class:`transfer <DarkAges.transfer.transfer>` based on the tables from XXXX.ZZZZ.
*  :code:`logEnergies`: Array containing the logarithm of the kinetic energies of
   the particles to the base 10 (The energies are given in units of
   :math:`\\mathrm{eV}`). Per default this are the same energies as in the
   :code:`transfer_functions`.
   This array can be accesed by the getter and setter functions :meth:`set_logEnergies`
   and :meth:`get_logEnergies`.
*  :code:`redshift`: Array containing the values of the redshifts (:math:`z+1`)
   at which the energy is *injected* into the primordial plasma.
   As of now this array equals the redshift-array in :code:`transfer_functions`,
   hence only the getter-function :code:`get_redshift` can be used.
"""

from __future__ import absolute_import, division, print_function
from builtins import range

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
transfer_functions_corr = None
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
		Masked message to be printed in :obj:`stdout`.
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
	def __init__(self, message, reason=None):
		Exception.__init__(self)
		self.message = message
		self.reason = reason

		self.name = self.__class__.__name__

	def __str__(self):
		if self.reason is not None:
			return '\n\n!!! ERROR ({} - Reason: {}) !!!\n\n --> {} \n'.format(self.name,self.reason,self.message)
		else:
			return '\n\n!!! ERROR ({}) !!!\n\n --> {} \n'.format(self.name,self.message)



def get_background(key=None):
	u"""Returns the current global values of the cosmological background
	parameters :math:`H_0` *(in 1/s)*, :math:`\\Omega_\\mathrm{matter}`
	and :math:`\\Omega_\\mathrm{radiation}`.

	The return value is either the complete :code:`CosmoBackground` dictionary

	.. code::

		CosmoBackground = {'H0': ...,'Omega_m': ...,'Omega_r':...}

	when no key is given or only one entry when the respective key (:code:`H0`,
	:code:`Omega_m` or :code:`Omega_r`) is given as an argument of the function.

	Parameters
	----------
	key : :obj:`str`, *optional*
		Key to return only a specific entry of the :code:`CosmoBackground` dictionary.
		If not given, the complete dictionarz is returned.

	Returns
	-------
	:obj:`float` or :obj:`dict`
		Value of :code:`CosmoBackground[key]`, if key is not None. Else Returns
		the complete dictionary.
	"""
	if CosmoBackground is None:
		raise DarkAgesError('The background parameters are not set.')
	if key is None:
		return CosmoBackground
	else:
		if key in CosmoBackground:
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
	u"""Returns the global array with the values of :math:`z+1` used at various
	places throughout the code.

	Returns
	-------
	:obj:`array-like`
		Array with values of the redshift (:math:`z+1`)
	"""

	if redshift is None:
		raise DarkAgesError('"redshift" was not initialized yet.')
	else:
		return redshift

def set_redshift(z):
	u"""Feed the global array :code:`redshift` with the values given in
	the input :code:`z`.

	.. note::
		Currently it is only sensible to use the same values for the redshifts
		as defined by the grid of the transfer functions. For that reason this
		function can only be run once, when the global array :code:`redshift`
		is not yet initialized.

	Parameters
	----------
	z : :obj:`array-like`
		Array of values of the redshift :math:`z+1` which should be stored as
		the globally accessible array :code:`redshift`.
	"""

	global redshift
	# As of now it makes only sense to use the redshift of the transfer functions
	# rather than a custom array, because the interpolation of the transfer functions
	# in redshift space should be done with caution and is not included in this
	# Version of the module. Therefore we only allow that this function can
	# be exectued once when the initial value of redshift is set.
	if redshift is None:
		redshift = z

def get_logEnergies():
	u"""Returns the global array with the values of :math:`\\log_{10} (E)` used
	at various places throughout the code.

	Returns
	-------
	:obj:`array-like`
		Array with values of :math:`\\log_{10} (E)`.
	"""

	if logEnergies is None:
		raise DarkAgesError('"logEnergies" was not initialized yet.')
	else:
		return logEnergies

def set_logEnergies(logE):
	u"""Feed the global array :code:`logEnergies` with the values given in
	the input :code:`logEnergies`.

	The input contains the logarithm of the energies at which the spectra and the
	transfer functions should be sampled. The base is 10 and the energies are given
	in units of :math:`\\mathrm{eV}`.

	Parameters
	----------
	logE : :obj:`array-like`
		Array of values of :math:`\\log_{10} (E)` which should be stored as
		the globally accessible array :code:`logEnergies`.
	"""

	global logEnergies
	logEnergies = logE

def _transfer_init_and_dump():
	global transfer_functions
	global transfer_functions_corr
	for channel in list(channel_dict.keys()):
		idx = channel_dict.get(channel)
		transfer_functions[idx] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch{:d}.dat'.format(idx+1)))
		transfer_dump(transfer_functions[idx], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch{:d}.obj'.format(idx+1)))
	transfer_functions_corr = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Corr.dat'))
	transfer_dump(transfer_functions_corr, os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Corr.obj'))

def _transfer_load_from_dump():
	global transfer_functions
	global transfer_functions_corr
	for channel in list(channel_dict.keys()):
		idx = channel_dict.get(channel)
		transfer_functions[idx] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch{:d}.obj'.format(idx+1)) )
	transfer_functions_corr = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Corr.obj') )

#################################

if (transfer_functions is None) or (transfer_functions_corr is None):
	transfer_functions = np.empty(shape=5, dtype=transfer)

	transfer_is_initialized = True
	for i in range(5):
		transfer_is_initialized = transfer_is_initialized and os.path.isfile(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch{:d}.obj'.format(i+1)))
	transfer_is_initialized = transfer_is_initialized and os.path.isfile(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Corr.obj'))

	if not transfer_is_initialized:
		print_info('The transfer seem not to be initialized. This will be done now. this may take a few seconds.')
		_transfer_init_and_dump()
		print_info('The transfer functions are now initialized and loaded.\n')
	else:
		#print_info('The transfer functions are already initialized and loaded.\n')
		_transfer_load_from_dump()
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
