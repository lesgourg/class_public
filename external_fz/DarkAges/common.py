u"""
.. module:: common
   :synopsis: Provide helper functions for the calculation of f(z)
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

Collection of functions needed to calculate the energy deposition.
"""

from scipy.integrate import trapz
from scipy.special import erf
import os
import sys
import numpy as np

###############################
### Definition of functions ###
###############################

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

def set_cosmology(H0 = 67.27, Om_M = 0.315, Om_R = 8e-5):
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

	Cosmo_H0 = H0 * _km_per_Mpc # in 1/s
	Cosmo_Omega_M = Om_M
	Cosmo_Omega_R = Om_R
	return Cosmo_H0, Cosmo_Omega_M, Cosmo_Omega_R

if ('DARKAGES_BG_H0' in os.environ) and ('DARKAGES_BG_OM_M' in os.environ) and ('DARKAGES_BG_OM_R' in os.environ):
	Cosmo_H0, Cosmo_Omega_M, Cosmo_Omega_R = set_cosmology(H0 = float(os.environ['DARKAGES_BG_H0']),
                                                           Om_M = float(os.environ['DARKAGES_BG_OM_M']),
                                                           Om_R = float(os.environ['DARKAGES_BG_OM_R']))
else:
	Cosmo_H0, Cosmo_Omega_M, Cosmo_Omega_R = set_cosmology()

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

def logConversion( log_array , base=10 ):
	u"""Returns :code:`pow(base, entry)` on every entry in :code:`log_array`

	Parameters
	----------
	log_array : :obj:`array-like`
		Array (:code:`shape = (k)`) of exponents to the base :code:`base`
	base : :obj:`float`, :obj:`int`, *optional*
	 	base to which the exponents in :code:`log_array` are defined to.
		If not given the base 10 is assumed

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k)`) on which :code:`pow(base, entry)` for
		every entry in :code:`log_array` was applied on.
	"""

	def dummy(single_log10_num):
		return base**single_log10_num
	return np.vectorize(dummy).__call__(log_array)

def simple_H(redshift, H0 = Cosmo_H0, Omega_M = Cosmo_Omega_M):
	u"""Returns the Hubble parameter at given redshift in a purely matter
	dominated universe.

	.. note::
		This function is not used in the standard calculations of this package

	Parameters
	----------
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) of values :math:`z+1`
	H0 : :obj:`float`, *optional*
		Todays Hubble parameter **(in units of 1/s)**. If not given the standard
		PLANCK-bestfit value (:math:`H_0 = 67.27\\; \\mathrm{km} /\\mathrm{Mpc\\,s}`) is assumed.
	Omega_M : :obj:`float`, *optional*
		Todays matter fraction. If not given the standard PLANCK-bestfit value
		(:math:`\\Omega_\\mathrm{matter} = 0.315`) is assumed.

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k)`) of Hubble parameters at the redshifts given in :code:`redshift`.
	"""

	return H(redshift, Omega_R = 0.)

def H(redshift, H0=Cosmo_H0, Omega_M = Cosmo_Omega_M, Omega_R = Cosmo_Omega_R):
	u"""Returns the Hubble parameter at given redshift

	Parameters
	----------
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) of values :math:`z+1`
	H0 : :obj:`float`, *optional*
		Todays Hubble parameter **(in units of 1/s)**. If not given the standard
		PLANCK-bestfit value (:math:`H_0 = 67.27\\; \\mathrm{km} /\\mathrm{Mpc\\,s}`) is assumed.
	Omega_M : :obj:`float`, *optional*
		Todays matter fraction. If not given the standard PLANCK-bestfit value
		(:math:`\\Omega_\\mathrm{matter} = 0.315`) is assumed.
	Omega_R : :obj:`float`, *optional*
		Todays radiation fraction. If not given the standard PLANCK-bestfit value
		(:math:`\\Omega_\\mathrm{radiation} = 8\\cdot10^{-5}`) is assumed.

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k)`) of Hubble parameters at the redshifts given in :code:`redshift`.
	"""

	return H0 * np.sqrt( redshift**(3.) * Omega_M + redshift**(4.) * Omega_R + (1-Omega_R-Omega_M) )
	#return H0 * np.sqrt( redshift**(3.) * Omega_M + redshift**(4.) * Omega_R ) * 1/np.sqrt(Omega_R + Omega_M)

def time_at_z(redshift, H0=Cosmo_H0, Omega_M = Cosmo_Omega_M, Omega_R = Cosmo_Omega_R):
	u"""Returns time (in seconds) at a given redshift.

	For simplicity it is assumed that only matter and radiation are present
	and dark energy is negligible. Valid for high redshifts
	(Taken from `EnergyAbsorptionCalculator.nb` provided as a supplement of
	`arXiV:1211.0283 <https://arxig.org/abs/1211.0283>`_)

	Parameters
	----------
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) of values :math:`z+1`
	H0 : :obj:`float`, *optional*
		Todays Hubble parameter (in units of 1/s). If not given the standard
		PLANCK-bestfit value (:math:`H_0 = 67.27\\; \\mathrm{km} /\\mathrm{Mpc\\,s}`) is assumed.
	Omega_M : :obj:`float`, *optional*
		Todays matter fraction. If not given the standard PLANCK-bestfit value
		(:math:`\\Omega_\\mathrm{matter} = 0.315`) is assumed.
	Omega_R : :obj:`float`, *optional*
		Todays radiation fraction. If not given the standard PLANCK-bestfit value
		(:math:`\\Omega_\\mathrm{radiation} = 8\\cdot10^{-5}`) is assumed.

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k)`) of t (in seconds) at the redshifts given in :code:`redshift`.
	"""

	return 2/(3 * Omega_M**2 * redshift * H0) * ( Omega_M * np.sqrt(Omega_R + (Omega_M / redshift)) + 2 * Omega_R**1.5 * redshift - 2 * Omega_R * np.sqrt(redshift*(Omega_M + redshift*Omega_R) ) )

def conversion( redshift, alpha=3 ):
	u"""Returns :math:`\\frac{\\left(z+1\\right)^\\alpha}{H(z)}`

	Returns conversion factor in front of the z and energy integrals
	which takes the expansion of comoving volume into account.

	The exponent :math:`alpha` depends on the particular energy injection
	history you consider.(3 for annihilating DM (s-waves), 0 for a decaying species)

	Parameters
	----------
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) of values :math:`z+1`
	alpha : :obj:`int`, :obj:`float`, *optional*
		Exponent to take the scaling of the number density in comoving volume
		into account (e.g. 3 for s-wave-annihilation, 0 for decay).
		If not given, s-wave-annihilation (:math:`\\alpha = 3`) is assumed.

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k)`) of the conversion factor at the redshifts given in :code:`redshift`.
	"""

	return ( redshift**alpha / H(redshift) )

def nan_clean( input_array ):
	u"""Returns the input where all invalid entries (NaN, infinity, ..) are
	replaced with zero.

	Parameters
	----------
	input_array : :obj:`array-like`
		Array with potentially invalid numerical entries
		(:code:`NaN`, :code:`infinity`,..)

	Returns
	-------
	:obj:`array-like`
		Copy of :code:`input_array` with all invalid entries replaced by zero.
	"""

	def dummy( single_input ):
		if abs(single_input) < np.inf:
			return single_input
		else:
			return 0.
	return np.vectorize(dummy).__call__(input_array)

def get_index( array, entry ):
	u"""Returns the index of the first occurence of a specific entry in a given array.

	Parameters
	----------
	array : :obj:`array-like`
		Array to search in
	entry : :obj:`float`, :obj:`int`
		Entry to search for

	Returns
	-------
	:obj:`int`
		Index of the first occurence of :code:`entry` in :code:`array`.
	"""
	return np.where(array == entry)[0][0].astype(np.int32)


def boost_factor_halos(redshift,zh,fh):
	ret = 1 + fh*erf(redshift/(1+zh))/redshift**3
	return ret

def scaling_boost_factor(redshift,spec_point,zh,fh):
	ret = spec_point*(1 + fh*erf(redshift/(1+zh))/redshift**3)
	return ret

def f_function(logE, z_inj, z_dep, normalization,
               transfer_phot, transfer_elec,
               spec_phot, spec_elec, alpha=3, **kwargs):
	u"""Returns the effective efficiency factor :math:`f_c (z)`
	for the deposition channel :math:`c`.

	This method calculates the effective efficiency factor in dependence
	of the redshift

	......WRITE HERE .....

	Parameters
	----------
	logE : :obj:`array-like`
		Array (:code:`shape = (l)`) of the logarithms of the kinetic energies of teh particles
		(in units of eV) to the base 10.
	z_inj : :obj:`array-like`
		Array (:code:`shape = (m)`) of the values :math:`z_\\mathrm{inj.}` at which the energy
		was injected (e.g. by annihilation or decay)
	z_dep : :obj:`array-like`
		Array (:code:`shape = (k)`) of the values :math:`z_\\mathrm{dep.}` at which the energy
		was deposited into the IGM
	normalization : :obj:`array-like`
		Array (:code:`shape = (k)`) containing the proper normalization of the injected spectra
		of photons and electrons at each timestep / at each redshift of deposition
	transfer_phot : :obj:`array-like`
		Array (:code:`shape = (k,l,m)`) containing the discretized transfer functions
		:math:`T^\\mathrm{phot.}_{klm}` for photons
	transfer_elec : :obj:`array-like`
		Array (:code:`shape = (k,l,m)`) containing the discretized transfer functions
		:math:`T^\\mathrm{elec.}_{klm}` for electrons and positrons
	spec_phot : :obj:`array-like`
		Array (:code:`shape = (l,m)`) containing the double differential spectrum
		:math:`\\frac{\\mathrm{d}^2 N}{ \\mathrm{d}E \\mathrm{d}t }` of photons
	spec_elec : :obj:`array-like`
		Array (:code:`shape = (l,m)`) containing the double differential spectrum
		:math:`\\frac{\\mathrm{d}^2 N}{ \\mathrm{d}E \\mathrm{d}t }` of electrons
		and positrons.
	alpha : :obj:`int`, :obj:`float`, *optional*
		Exponent to take the scaling of the number density in comoving volume
		into account (see also: :meth:`conversion <DarkAges.common.conversion>`).
		If not given the default value :math:`\\alpha = 3` is taken.

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k)`) of :math:`f_c (z)` at the redshifts of
		deposition given in :code:`z_dep`
	"""

	E = logConversion(logE)

	how_to_integrate = kwargs.get('E_integration_scheme','logE')
	if how_to_integrate not in ['logE','energy']:
		print_error('The energy integration-scheme >> {0} << is not known'.format(how_to_integrate))

	norm = ( conversion(z_dep,alpha=alpha) )*( normalization )

	energy_integral = np.zeros( shape=(len(z_dep),len(z_inj)), dtype=np.float64)
	for i in xrange(len(energy_integral)):
		if how_to_integrate == 'logE':
			for k in xrange(i,len(energy_integral[i])):
				int_phot = transfer_phot[i,:,k]*spec_phot[:,k]*(E**2)/np.log10(np.e)
				int_elec = transfer_elec[i,:,k]*spec_elec[:,k]*(E**2)/np.log10(np.e)
				energy_integral[i][k] = trapz( int_phot + int_elec, logE )
		elif how_to_integrate == 'energy':
			for k in xrange(i,len(energy_integral[i])):
				int_phot = transfer_phot[i,:,k]*spec_phot[:,k]*(E**1)
				int_elec = transfer_elec[i,:,k]*spec_elec[:,k]*(E**1)
				energy_integral[i][k] = trapz( int_phot + int_elec, E )

	z_integral = np.zeros_like( z_dep, dtype=np.float64)
	dummy = np.arange(1,len(z_inj)+1)
	for i in xrange(len(z_integral)):
		low = max(i,0)
		#low = i
		integrand = ( conversion(z_inj[low:], alpha=alpha) )*energy_integral[i,low:]
		z_integral[i] = trapz( integrand, dummy[low:] )

	result = np.empty_like( norm, dtype=np.float64 )
	for i in xrange(len(norm)):
		if norm[i] != 0 and abs(z_integral[i]) < np.inf :
			result[i] = (z_integral[i] / norm[i])
		else:
			#result[i] = np.nan
			result[i] = 0

	return result

def log_fit(points,func,xgrid,exponent=1):
	u"""Returns an array of interpolated points using the
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`-class.

	This method is used to initialize and call the
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`
	at the same time. It is used if the interpolated function needs to be
	read imediately rather than to be stored for later use.

	Parameters
	----------
	points : :obj:`array-like`
		Array (:code:`shape = (k)`) with the points at which the fuction to
		interpolate is given.
	func : :obj:`array-like`
		Array (:code:`shape = (k)`) with the values of the function to interpolate.
	xgrid : :obj:`array-like`
		Array (:code:`shape = (l)`) with the points at which teh function should
		be interpolated. Needs to fulfill :code:`min(xgrid) >= min(points)` and
		:code:`max(xgrid) <= max(points)`.
	exponent : :obj:`int`, :obj:`float`, *optional*
		Exponent to specify the powers of :code:`points` mulitplied to
		:code:`func` before the function is transformed into logspace.
		(see also: :class:`logInterpolator <DarkAges.interpolator.logInterpolator>`).
		If not given, :code:`points` is multiplied linearly (:code:`exponent=1`).

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (l)`) with the interpolated values of the function
		at the points given by :code:`xgrid`.
	"""

	from .interpolator import logInterpolator
	tmp_interpolator = logInterpolator(points, func, exponent)
	#tmp_interpolator = logLinearInterpolator(points, func, exponent)
	out = tmp_interpolator(xgrid)
	return out

def sample_spectrum(input_spec_el, input_spec_ph, input_spec_oth, input_log10E, m, sampling_log10E, **kwargs):
	u"""Returns the interpolated and properly normalized particle spectrum

	This method interpolates the particle spectra defined at the points
	:code:`input_log10E`, applies the normalization given the injection history
	in question and returns the recurrent spectra ath the points given in
	:code:`sampling_log10E`

	This is mainly used to interpolate the reference spectra for the use
	of the models :class:`annihilating_model <DarkAges.model.annihilating_model>`
	or :class:`decaying_model <DarkAges.model.decaying_model>`.
	In the most usual case the reference spectra are read from a table
	with the kinetic energies (i.e. their logarithm to the base 10) given in
	:code:`input_log10E` which do not coincide with the energies of the
	:class:`transfer_instance <DarkAges.transfer.transfer>` given in
	:code:`sampling_log10E`.
	In this method the interpolation is done using DarkAges's own
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`
	(with :code:`exponent = 1`).
	Any unphysical interpolated point (e.g negative value of the function or
	non-zero value for a kinetic enrgy of the particle higher than the rest mass
	of the initial DM candidate) will be set to zero.

	Parameters
	----------
	input_spec_el : :obj:`array-like`
		Array (:code:`shape = (k)`) of the differential spectrum
		:math:`\\frac{\\mathrm{d}N}{\\mathrm{d}E}` of electrons and positrons.
	input_spec_ph : :obj:`array-like`
		Array (:code:`shape = (k)`) of the differential spectrum
		:math:`\\frac{\\mathrm{d}N}{\\mathrm{d}E}` of photons.
	input_spec_oth : :obj:`array-like`
		Array (:code:`shape = (k)`) of the differential spectrum
		:math:`\\frac{\\mathrm{d}N}{\\mathrm{d}E}` of particles not injecting
		energy into the IGM. Used for the proper normailzation of the spectra.
	input_log10E : :obj:`array-like`
		Array (:code:`shape = (k)`) of the logarithm of the kinetic energies
		of the particles to the base 10 at which the input spectra are
		defined.
	m : :obj:`float`
		Masss of the DM candidate.
	sampling_log10E : :obj:`array-like`
		Array (:code:`shape = (l)`) of the logarithm of the kinetic energies
		(*in units of* :math:`\\mathrm{eV}`) of the particles to the base 10
		at which the spectra should be interpolated.

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (3,l)`) of the sampled particle spectra of
		positrons and electrons, photons, and particle not interacting
		with the IGM at the energies specified in :code:`sampling_log10E`.
	"""

	scale = kwargs.get('scale','GeV')
	spec_type = kwargs.get('spec_type','dN/dE')
	hist = kwargs.get('injection_history','annihilation')
	norm_dict = {'annihilation':2., 'decay':1.}
	norm = norm_dict.get(hist,2.)*m

	failed = False
	try:
		assert spec_type in ['dN/dE','E.dN/dE']
	except AssertionError:
		print_warning('Unknown type of input spectrum >> {0} <<. Valid options are: "dN/dE" and "E.dN/dE"'.foramt(spec_type))
		failed = True

	try:
		assert scale in ['ev','keV','MeV','GeV']
	except AssertionError:
		print_warning('Unknown scale of your input >> {0} <<. Valid options are: "eV", "keV", "MeV", and "GeV"'.format(scale) )
		failed = True
		scale = 'eV' # Set an abitrary, but valid scale. Since the method is labeled as failed it will return zero anywway

	scale_dict = {
          'GeV': [1e9, 1e-9, 9.],
          'MeV': [1e3, 1e-3, 6.],
          'keV': [1e3, 1e-3, 3.],
           'eV': [1e0, 1e-0, 0.],
	}

	m = scale_dict[scale][0]*m
	norm *= scale_dict[scale][0]
	input_log10E += scale_dict[scale][2]*np.ones_like(input_log10E).astype(np.float64)
	# Check if the spectrum is normalized to "integral(E*dN/dE dE) = norm , where norm is chosen depending on the injection history"
	#    If not and the integral is non_zero, rescale the spectrum, if the spectrum is zero (e.g. being in a kinematically forbidden regioern of param space)
	#    return a zero_spectrum
	if spec_type == 'dN/dE':
		factor1 = logConversion(input_log10E)
		factor2 = np.ones_like(input_log10E).astype(np.float64)
	else:
		factor1 = np.ones_like(input_log10E).astype(np.float64)
		factor2 = 1. / logConversion(input_log10E)
	total_dep_energy = trapz( factor1*(input_spec_el+input_spec_ph+input_spec_oth)*scale_dict[scale][1], logConversion(input_log10E) )

	non_zero_spec = ( total_dep_energy > 0.)
	if non_zero_spec and not failed:
		rescaling = total_dep_energy / norm

		out_el = log_fit(input_log10E, factor2*scale_dict[scale][1]*input_spec_el/rescaling, sampling_log10E)
		out_ph = log_fit(input_log10E, factor2*scale_dict[scale][1]*input_spec_ph/rescaling, sampling_log10E)
		out_oth = log_fit(input_log10E, factor2*scale_dict[scale][1]*input_spec_oth/rescaling, sampling_log10E)

		# In the unphysical region E > m the spectrum needs to vanish
		unphysical_region_mask =  (logConversion(sampling_log10E) > m)
		out_el[unphysical_region_mask] = 0.
		out_ph[unphysical_region_mask] = 0.
		out_oth[unphysical_region_mask] = 0.
	else:
		out_el = np.zeros_like(sampling_log10E).astype(np.float64)
		out_ph = np.zeros_like(sampling_log10E).astype(np.float64)
		out_oth = np.zeros_like(sampling_log10E).astype(np.float64)

	return np.array([out_el, out_ph, out_oth])

def finalize(redshift, f_heat, f_lya, f_ionH, f_ionHe, f_lowE, **kwargs):
	u"""Prints the table of redshift and :math:`f_c(z)` for the deposition
	channels in question into :obj:`stdout`

	Since `CLASS <http://class-code.net>`_ expects a pure table of :math:`z` and
	:math:`f_c()z` this is the last function called during a DarkAges-session.

	.. warning::
		For the correct usage of this package together with
		`CLASS <http://class-code.net>`_ the only allowed output
		are line with a single number, containing the number of the lines
		of the table to follow and the table

		+-------+------+-----+-------+--------+------+
		|   #z  | Heat | LyA | H-Ion | He-Ion | LowE |
		+=======+======+=====+=======+========+======+
		|   0   |      |     |       |        |      |
		+-------+------+-----+-------+--------+------+
		|  ...  |  ... | ... |  ...  |   ...  |  ... |
		+-------+------+-----+-------+--------+------+
		| 10000 |      |     |       |        |      |
		+-------+------+-----+-------+--------+------+

		Please make sure that all other message printed are silenced or
		at least covered by '#' (see :meth:`print_info`)

	Parameters
	----------
	redshift : :obj:`array-like`
		Array (:code:`shape = (k)`) with the values of redshift :math:`z`.
		Note that here *redshift* is meant to be :math:`z` and not
		:math:`z+1`
	f_heat : :obj:`array-like`
		Array (:code:`shape = (k)`) with the values of the effective efficiency factor
		for the deposition channel *Haeting* at the redshifts given by :code:`redshift`
	f_lya : :obj:`array-like`
		As :code:`f_heat` but for the deposition channel *Ly-a excitation*
	f_ionH : :obj:`array-like`
		As :code:`f_heat` but for the deposition channel *hydrogen ionization*
	f_ionHe : :obj:`array-like`
		As :code:`f_heat` but for the deposition channel *helium ionization*
	f_lowE : :obj:`array-like`
		As :code:`f_heat` but for the deposition channel *sub10.2ev - photons*
	"""

	redshift = redshift - np.ones_like(redshift) # Go from DarkAges-redshift (z+1) to CLASS-redshift (z)

	first = int(kwargs.get('first_index',1))
	last_idx = int(kwargs.get('last_index',1))
	if last_idx == 0:
		print_warning('We strongly discourage you to assign the value "0" to "last_index". The last entry of the table is zero for numerical reasons.')
	last = len(redshift) - last_idx
	min_z = kwargs.get('lower_z_bound',0.)
	max_z = kwargs.get('upper_z_bound',1e4)
	print(50*'#')
	print('### This is the standardized output to be read by CLASS.\n### For the correct usage ensure that all other\n### "print(...)"-commands in your script are silenced.')
	print(50*'#'+'\n')
	print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n{:d}\n'.format( (last-first) + 2))
	print('{:.2e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}'.format(min_z,f_heat[first],f_lya[first],f_ionH[first],f_ionHe[first],f_lowE[first]))
	for idx in range(first,last):
		print('{:.2e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}'.format(redshift[idx],f_heat[idx],f_lya[idx],f_ionH[idx],f_ionHe[idx],f_lowE[idx]))
	print('{:.2e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}'.format(max_z,f_heat[last-1],f_lya[last-1],f_ionH[last-1],f_ionHe[last-1],f_lowE[last-1]))
