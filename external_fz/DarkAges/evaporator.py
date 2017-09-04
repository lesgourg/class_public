u"""
.. module:: evaporator
   :synopsis: Definition of functions in the scope of evaporating primordial black holes.
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

This module contains the definition of functions which are useful for the calculation
of the evaporation of an evaporating black holes. This are mainly the following
functions:

#. Functions related to the evolution of the mass of the primordial black holes
   with time:

	* :meth:`PBH_dMdt`
	* :meth:`PBH_mass_at_z` and :meth:`PBH_mass_at_t`

#. Functions related to the spectrum of the particles produced by the evaporation:

	* :meth:`PBH_spectrum_at_m`, :meth:`Primary_spectrum` and :meth:`full_spectrum`
	* :meth:`PBH_F_of_M`, :meth:`fraction_at_M`

#. Functions to translate between the mass of a black hole and its temperature

	* :meth:`get_mass_from_temperature` and :meth:`get_temperature_from_mass`

"""

import numpy as np

from .common import time_at_z, print_warning, logConversion, H, nan_clean
from .__init__ import redshift, logEnergies, DarkAgesError
from scipy.integrate import odeint as solve_ode
from scipy.integrate import trapz

#particle_dict has structure 'particle_name':[mass,spin_label,spin_orientations, flavorfactor, antifactor, colorfactor, sigmoid_factor]
particle_dict = {'gamma': [0,'1',2,1,1,1,0],
			  'neutrino': [0,'1/2 N',1,3,2,1,0],
			  'electron': [5.110e-6,'1/2 C',2,1,2,1,0],
				  'muon': [1.057e-1,'1/2 C',2,1,2,1,0],
				   'tau': [1.777,'1/2 C',2,1,2,1,0],
					'up': [2.2e-3,'1/2 C',2,1,2,3,1],
				  'down': [4.7e-3,'1/2 C',2,1,2,3,1],
				 'charm': [1.28,'1/2 C',2,1,2,3,1],
			   'strange': [9.6e-2,'1/2 C',2,1,2,3,1],
				   'top': [173.1,'1/2 C',2,1,2,3,1],
				'bottom': [4.18,'1/2 C',2,1,2,3,1],
			   'w-boson': [80.39,'1',3,1,2,1,0],
			   'z-boson': [91.19,'1',3,1,1,1,0],
				 'gluon': [6e-1,'1',2,8,1,1,1],
				 'higgs': [125.09,'0',1,1,1,1,0],
				   'pi0': [1.350e-1,'0',1,1,1,1,-1],
				  'piCh': [1.396e-1,'0',1,1,2,1,-1]
}

def _particle_list_resolver( *particles ):
	particles = list(particles)
	if ('ALL' in particles) or ('all' in particles):
		return particle_dict.keys()
	else:
		if ('light quarks' in particles):
			particles.pop(particles.index('light quarks'))
			particles.extend(['up','down','strange'])
		if ('pions' in particles):
			particles.pop(particles.index('pions'))
			particles.extend(['pi0','piCh'])
		return particles

def _get_spin(particle):
	spin_label = particle_dict.get(particle)[1]
	tmp_dict = {'0':0.0,
			'1/2 N':0.5,
			'1/2 C':0.5,
				'1':1.0}
	return tmp_dict.get(spin_label)

def PBH_spectrum_at_m( mass, logEnergies, *particles):
	ret = np.zeros((len(logEnergies),len(mass)), dtype=np.float64)
	E = logConversion(logEnergies - 9)

	if particles:
		particles = _particle_list_resolver( *particles )
		for particle in particles:
			spin = _get_spin(particle)
			fraction = fraction_at_M( mass, particle )
			ret[:,:] += fraction*np.asarray(np.vectorize(Primary_spectrum, excluded = ['spin']).__call__( energy=E[:,None], PBH_mass=mass[None,:], spin=spin))
	else:
		raise DarkAgesError('There is no particles given')
	return ret

#def PBH_spectrum( PBH_mass_ini, *particles):
#	ret = np.zeros((len(logEnergies), len(redshift)), dtype=np.float64)
#	mass_at_z = PBH_mass_at_z(PBH_mass_ini, z_start =1e4)[-1,:]
#	for idx, mass in enumerate(mass_at_z):
#		ret[:,idx,None] = PBH_spectrum_at_m( mass, logEnergies, *particles )
#	return ret

#def full_spectrum( PBH_mass_ini ):
#	ret = np.zeros((3, len(logEnergies), len(redshift)), dtype=np.float64)
#	mass_at_z = PBH_mass_at_z(PBH_mass_ini, z_start =1e4)[-1,:]
#	E = logConversion(logEnergies - 9)
#	ret[0,:,:] = PBH_spectrum_at_m( mass_at_z, logEnergies, 'electron')
#	ret[1,:,:] = PBH_spectrum_at_m( mass_at_z, logEnergies, 'gamma')
#	ret[2,:,:] = PBH_spectrum_at_m( mass_at_z, logEnergies, 'ALL')
#	return E, redshift, ret

def get_temperature_from_mass(mass):
	u"""Returns the temperature of a black hole of a given mass

	Parameters
	----------
	mass : :obj:`float`
		Mass of the black hole (*in units of* :math:`\\mathrm{g}`)

	Returns
	-------
	:obj:`float`
		Temperature of the black hole (*in units of* :math:`\\mathrm{GeV}`)
	"""

	return 1.06e13 / mass

def get_mass_from_temperature(temperature):
	u"""Returns the mass of a black hole of a given temperature

	Parameters
	----------
	temperature : :obj:`float`
		Temperature of the black hole (*in units of* :math:`\\mathrm{GeV}`)

	Returns
	-------
	:obj:`float`
		Mass of the black hole (*in units of* :math:`\\mathrm{g}`)
	"""

	return 1.06e13 / temperature

def Primary_spectrum( energy, PBH_mass, spin ):
	u"""Returns the double differential spectrum
	:math:`\\frac{\\mathrm{d}^2 N}{\\mathrm{d}E \\mathrm{d}t}` of particles with
	a given :code:`spin` and kinetic :code:`energy` for an evaporating black hole
	of mass :code:`PBH_mass`.

	Parameters
	----------
	energy : :obj:`float`
		Kinetic energy of the particle (*in units of* :math:`\\mathrm{GeV}`)
		produced by the evaporating black hole.
	PBH_mass : :obj:`float`
		Current mass of the black hole (*in units of* :math:`\\mathrm{g}`)
	spin : :obj:`float`
		Spin of the particle produced by the evaporating black hole (Needs
		to be a multiple of :math:`\\frac{1}{2}`, i.e. :code:`2 * spin` is assumed
		to have integer value)

	Returns
	-------
	:obj:`float`
		Value of :math:`\\frac{\\mathrm{d}^2 N}{\\mathrm{d}E \\mathrm{d}t}`
	"""

	if abs(PBH_mass) < np.inf and PBH_mass > 0.:
		PBH_temperature = get_temperature_from_mass(PBH_mass)
		Gamma = 27 * ( 6.70861e-39)**2 * (energy * (PBH_mass * 5.61e23))**2
		return (1/(2*np.pi)) * Gamma / (np.exp( energy / PBH_temperature ) - (-1)**int(np.round(2*spin)))
	else:
		return 0.

def fraction_at_M( PBH_mass, channel, *all_channels):
	u"""Returns the relative contribution of a given particle (:code:`channel`)
	being emmited through evaporation of a black hole with mass :code:`PBH_mass`

	If not given via :code:`*all_channels`, the relative contribution is calculated
	with respect to all channels.

	Parameters
	----------
	PBH_mass : :obj:`float`
		Current mass of the black hole (*in units of* :math:`\\mathrm{g}`)
	channel : :obj:`str`
		Name of the channel in question (For the names of the possible channels
		see 'here <DarkAges.evaporator.particle_dict>')
	*all_channels : list of :obj:`str`, *optional*
		List of the channels to which the relative contribution should be calculated.
		If not given, all possible channels are considered

	Returns
	-------
	:obj:`float`
		Fraction of particles in :code:`channel` being emmited via black hole Evaporation
	"""

	if all_channels:
		return  PBH_F_of_M( PBH_mass, channel ) / PBH_F_of_M( PBH_mass, *all_channels )
	else:
		return  PBH_F_of_M( PBH_mass, channel ) / PBH_F_of_M( PBH_mass )

def PBH_F_of_M( PBH_mass, *particles ):
	u"""Returns the *effective degrees of freedom* :math:`F(M)` of a given list of
	particles of an evaporating black hole with mass :code:`PBH_mass`.

	If the list in :code:`*particles` contains the string :code:`'ALL'` then
	all particles are considered

	Parameters
	----------
	PBH_mass : :obj:`float`
		Current mass of the black hole (*in units of* :math:`\\mathrm{g}`)
	*particles : :obj:`str` or list of :obj:`str`, *optional*
		Particles to include in the calculation of :math:`\\matrhm{F}(M)`. If
		left blank, all particles are considered.
	Returns
	-------
	:obj:`float`
		Value of :math:`F(M)`
	"""

	# factor_dict has the structure 'spin_label' : [f, beta]
	factor_dict = {'0':[0.267, 2.66],
			   '1/2 N':[0.142, 4.53],
			   '1/2 C':[0.147, 4.53],
			       '1':[0.060, 6.04]}

	def _single_contribution(PBH_mass, mass_of_particle, multiplicity, beta, simoid_factor = 0):
		if sigmoid_factor != 0:
			T_BH = get_temperature_from_mass(PBH_mass)
			Lam_QH = 0.3
			sigma = 0.1
			activation = 1. / (1. + np.exp(- sigmoid_factor*(np.log10(T_BH) - np.log10(Lam_QH))/sigma))
			#activation = np.heaviside(sigmoid_factor*(T_BH - Lam_QH), 0.5)
		else:
			activation = 1
		if mass_of_particle > 0:
			equivalent_BH_mass = get_mass_from_temperature(mass_of_particle)
			return activation * multiplicity * np.exp(-PBH_mass / (beta*equivalent_BH_mass) )
		else:
			return activation * multiplicity

	ret = 0.
	if not particles:
		particles = _particle_list_resolver( 'ALL' )
	else:
		particles = _particle_list_resolver( *particles )
	for particle in particles:
		if particle in particle_dict:
			values = particle_dict.get(particle)
			mass_of_particle = values[0]
			multiplicity = factor_dict[values[1]][0] * values[2] * values[3] * values[4] * values[5]
			beta = factor_dict[values[1]][1]
			sigmoid_factor = values[6]
			ret += _single_contribution(PBH_mass, mass_of_particle, multiplicity, beta, sigmoid_factor)
		else:
			print_warning('The particle "{:s}" is not recognized'.format(particle))

	return ret

def PBH_dMdt(PBH_mass, time, scale=1):
	u"""Returns the differential mass loss rate of an primordial black hole with
	a mass of :code:`PBH_mass*scale` gramm.

	This method calculates the differential mass loss rate at a given (rescaled) mass
	of the black hole and at a given time, given by

	.. math::

	   \\frac{\\mathrm{d}M\,\\mathrm{[g]}}{\\mathrm{d}t} = - 5.34\\cdot 10^{25}
	   \\cdot \\mathcal{F}(M) \\cdot \\left(\\frac{1}{M\,\\mathrm{[g]}}\\right)^2
	   \,\\frac{\\mathrm{1}}{\\mathrm{s}}

	.. note::

	   Even if this method is not using the variable :code:`time`, it is needed for the
	   ODE-solver within the :meth:`PBH_mass_at_t <DarkAges.evaporator.PBH_mass_at_t>`
	   for the correct interpretation of the differential equation for the mass of
	   the black hole.

	Parameters
	----------
	PBH_mass : :obj:`float`
		(Rescaled) mass of the black hole in units of :math:`\\mathrm{scale} \\cdot \\mathrm{g}`
	time : :obj:`float`
		Time in units of seconds. (Not used, but needed for the use with an ODE-solver)
	scale : :obj:`float`, *optional*
		For a better numerical performance the differential equation can be expressed in
		terms of a different scale than :math:`\\mathrm{g}`. For example the choice
		:code:`scale = 1e10` returns the differential mass loss rate in units of
		:math:`10^{10}\\mathrm{g}`. This parameter is optional. If not given, :code:`scale = 1`
		is assumed.

	Returns
	-------
	:obj:`float`
		Differential mass loss rate in units of :math:`\\frac{\\mathrm{scale}*\\mathrm{g}}{\\mathrm{s}}`
	"""

	if PBH_mass > 0:
		ret = -5.34e25*(scale**(-3))*PBH_F_of_M( scale*PBH_mass ) * (PBH_mass)**(-2)
		return ret
	else:
		return 0

#def PBH_mass_at_z(initial_PBH_mass, t_start = 1e10):
#	u"""Solves the differential equation for the black hole mass with the starting
#	conditions of an initial mass of :code:`initial_PBH_mass` at the initial
#	time of :code:`t_start`
#	"""

def PBH_mass_at_z(initial_PBH_mass, redshift=None, z_start = 1e4):
	def jac(mass, time, scale):
		out = np.zeros((1,1))
		out[0,0] = 0 # partial_(dMdt) / partial_t
		return out

	if redshift is None:
		from .__init__ import redshift as ham
		redshift = ham

	if np.max(redshift) >= z_start:
		raise DarkAgesError("The redshift array is in conflict with the redshift at which to start the evolution of the balck hole. At least one entry in 'redshift' exceeds the value of 'z_start'")

	log10_ini_mass = np.log10(initial_PBH_mass)
	scale = 10**(np.floor(log10_ini_mass)+5)
	initial_PBH_mass *= 1/scale

	temp_t = 10**np.linspace(np.log10(time_at_z(z_start)), np.log10(time_at_z(1.)), 1e5)
	temp_mass, full_info = solve_ode(PBH_dMdt, initial_PBH_mass, temp_t, Dfun=jac, args=(scale,), full_output=1,mxstep=10000)

	times =  time_at_z(redshift)
	PBH_mass_at_t = np.interp(times, temp_t, temp_mass[:,0])

	out = np.array([redshift, scale*PBH_mass_at_t])
	for idx in reversed(xrange(out.shape[-1])):
		if out[1,idx] <= 0:
			out[1,:idx+1] = 0.
			break

	return out
