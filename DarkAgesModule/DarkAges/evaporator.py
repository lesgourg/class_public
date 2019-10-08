u"""
.. module:: evaporator
   :synopsis: Definition of functions in the scope of evaporating primordial black holes.
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

This module contains the definition of functions which are useful for the calculation
of the evaporation of an evaporating black holes. This are mainly the following
functions:

#. Functions to translate between the mass of a black hole and its temperature

	* :meth:`get_mass_from_temperature <DarkAges.evaporator.get_mass_from_temperature>`
	  and :meth:`get_temperature_from_mass <DarkAges.evaporator.get_temperature_from_mass>`

#. Functions related to the spectrum of the particles produced by the evaporation:

	* :meth:`PBH_spectrum_at_m <DarkAges.evaporator.PBH_spectrum_at_m>`
	  and :meth:`PBH_primary_spectrum <DarkAges.evaporator.PBH_primary_spectrum>`
	* :meth:`PBH_F_of_M <DarkAges.evaporator.PBH_F_of_M>`
	  and :meth:`PBH_fraction_at_M <DarkAges.evaporator.PBH_fraction_at_M>`

#. Functions related to the evolution of the mass of the primordial black holes
   with time:

	* :meth:`PBH_dMdt <DarkAges.evaporator.PBH_dMdt>`
	* :meth:`PBH_mass_at_z <DarkAges.evaporator.PBH_mass_at_z>`

"""

from __future__ import absolute_import, division, print_function

import numpy as np

from .common import time_at_z, logConversion, H, nan_clean
from .__init__ import DarkAgesError, get_redshift
from scipy.integrate import odeint as solve_ode
from scipy.integrate import trapz

#particle_dict has structure 'particle_name':[mass,spin_label,spin_orientations, flavorfactor, antifactor, colorfactor, sigmoid_factor]
particle_dict = {'gamma': [0,'1',2,1,1,1,0],
			  'neutrino': [0,'1/2 N',1,3,2,1,0],
			  'electron': [5.110e-4,'1/2 C',2,1,2,1,0],
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
		return list(particle_dict.keys())
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

def PBH_F_of_M( PBH_mass, *particles, **DarkOptions ):
	u"""Returns the *effective degrees of freedom* :math:`F(M)` of a given list of
	particles of an evaporating black hole with mass :code:`PBH_mass`.

	If the list in :code:`*particles` contains the string :code:`'ALL'` then
	all particles are considered

	Parameters
	----------
	PBH_mass : :obj:`float`
		Current mass of the black hole (*in units of* :math:`\\mathrm{g}`)
	*particles : :obj:`str` or list of :obj:`str`, *optional*
		Particles to include in the calculation of :math:`\\mathcal{F}(M)`. If
		left blank, all particles are considered.
	Returns
	-------
	:obj:`float`
		Value of :math:`\\mathcal{F}(M)`
	"""

	# factor_dict has the structure 'spin_label' : [f, beta]
	factor_dict = {'0':[0.267, 2.66],
			   '1/2 N':[0.147, 4.53],
			   '1/2 C':[0.142, 4.53],
			       '1':[0.060, 6.04]}

	def _single_contribution(PBH_mass, mass_of_particle, multiplicity, beta, sigmoid_factor = 0, **DarkOptions):
		if sigmoid_factor != 0:
			use_QCD_transition = DarkOptions.get('QCD_phase_transition',True)
			T_BH = get_temperature_from_mass(PBH_mass)
			Lam_QCD = DarkOptions.get('QCD_lambda',0.3)
			if Lam_QCD <= 0:
				raise DarkAgesError('"QCD_lambda" needs to be positive, but yout input was >> {:.3e} <<'.format(Lam_QCD))
			sigma = DarkOptions.get('QCD_width',0.1)
			if sigma > 0.:
				if use_QCD_transition:
					activation = 1. / (1. + np.exp(- sigmoid_factor*(np.log10(T_BH) - np.log10(Lam_QCD))/sigma))
				elif sigmoid_factor < 0:
					activation = 0.
				else:
					activation = 1.
			else:
				raise DarkAgesError('Your input for "QCD_width" >> {:.3e} << is not valid. It needs to be a positive number'.format(sigma))
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
			ret += _single_contribution(PBH_mass, mass_of_particle, multiplicity, beta, sigmoid_factor, **DarkOptions)
		else:
			raise DarkAgesError('The particle "{:s}" is not recognized'.format(particle))

	return ret

def PBH_primary_spectrum( energy, PBH_mass, spin, **DarkOptions ):
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

	threshold_factor = DarkOptions.get('PBH_spec_threshold',3.)

	if abs(PBH_mass) < np.inf and PBH_mass > 0. and energy > 0.:
		PBH_temperature = get_temperature_from_mass(PBH_mass)
		Gamma = 27 * ( 6.70861e-39)**2 * (energy * (PBH_mass * 5.61e23))**2
		if energy >= threshold_factor*PBH_temperature:
			return (1/(2*np.pi)) * Gamma / (np.exp( energy / PBH_temperature ) - (-1)**int(np.round(2*spin)))
		else:
			return 0.
	else:
		return 0.

def PBH_fraction_at_M( PBH_mass, channel, *all_channels, **DarkOptions):
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
		return  PBH_F_of_M( PBH_mass, channel, **DarkOptions ) / PBH_F_of_M( PBH_mass, *all_channels, **DarkOptions )
	else:
		return  PBH_F_of_M( PBH_mass, channel, **DarkOptions ) / PBH_F_of_M( PBH_mass, **DarkOptions )

def PBH_spectrum_at_m( mass, logEnergies, *particles, **DarkOptions):
	u"""Returns the (combined) spectrum :math:`\\frac{\\mathrm{d}N}{\\mathrm{d}E}` of
	the particles given by the list :code:`*particles` emmited by a
	a black hole of mass :code:`mass` with a kinetic energy given by
	:code:`10**logEnergies`.

	This function computes for every particle the primary spectrum by
	:meth:`PBH_primary_spectrum <DarkAges.evaporator.PBH_primary_spectrum>` and the
	realtive contribution to the degrees of freedom :math:`\\mathcal{F}(M)` by
	:meth:`PBH_fraction_at_M <DarkAges.evaporator.PBH_fraction_at_M>` and computes the
	total spectrum.

	Parameters
	----------
	mass : :obj:`array-like`
		Array (:code:`shape = (k)`) of the black hole masses (*in units of* :math:`\\mathrm{g}`)
	logEnergies : :obj:`array-like`
		Array (:code:`shape = (l)`) of the logarithms (to the base 10) of the
		kinetic energy of the particles in question (the energy is given in units of GeV).
	*particles : tuple of :obj:`str`
		List of particles which should be considered. The contributions for each particle
		will be summed up. **Must at least contain one entry**

	Raises
	------
	DarkAgesError
		If no particles are given, (:code:`*particles` is empty)

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (k,l)`) of the summed particle spectrum (*in units of* :math:`\\mathrm{GeV}^{-1}`)
		at the enrgies and masses given by the inputs.
	"""

	ret = np.zeros((len(logEnergies),len(mass)), dtype=np.float64)
	E = logConversion(logEnergies - 9)

	if particles is not None:
		particles = _particle_list_resolver( *particles )
		for particle in particles:
			spin = _get_spin(particle)
			fraction = PBH_fraction_at_M( mass, particle, **DarkOptions )
			#energy = E
			energy = E + particle_dict.get(particle)[0]*np.ones_like(E)
			ret[:,:] += fraction*np.asarray(np.vectorize(PBH_primary_spectrum, excluded = ['spin']).__call__( energy=energy[:,None], PBH_mass=mass[None,:], spin=spin, **DarkOptions))
	else:
		raise DarkAgesError('There is no particle given')
	return ret

def PBH_dMdt(PBH_mass, time, scale=1, **DarkOptions):
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
		ret = -5.34e25*(scale**(-3))*PBH_F_of_M( scale*PBH_mass, **DarkOptions ) * (PBH_mass)**(-2)
		return ret
	else:
		return -0.0

def PBH_mass_at_z(initial_PBH_mass, redshift=None, **DarkOptions):
	u"""Solves the ODE for the PBH mass (:meth:`PBH_dMdt <DarkAges.evaporator.PBH_dMdt>`)
	and returns the masses at the redshifts given by the input :code:`redshift`

	If not specified by an additional keyword-argument in :code:`**DarkOptions`
	(by :code:`Start_evolution_at = ...`) the evolution is started at a redshift of 10.000

	Parameters
	----------
	initial_PBH_mass : :obj:`float`
		Initial mass of the primordial black hole (*in units of g*)
	redshift : :obj:`array-like` *optional*
		Array (:code:`shape = (l)`) of redshifts at which the PBH mass should
		be returned. If not given, the global redshift-array from
		:class:`the initializer <DarkAges.__init__>` is taken

	Returns
	-------
	:obj:`array-like`
		Array (:code:`shape = (l)`) of the PBH mass at the redshifts given in t
	"""

	# Jacobian of the ODE for the PBH mass.
	# Needed for better performance of the ODE-solver.
	def jac(mass, time, scale=1, **DarkOptions):
		out = np.ones((1,1))*(-2./(scale*mass))*PBH_dMdt(mass, time, scale=scale, **DarkOptions) # partial_(dMdt) / partial_M = -2*(dMdt)/M.
		#out = np.zeros((1,1))
		return out

	if redshift is None:
		redshift = get_redshift()

	z_start = DarkOptions.get('Start_evolution_at', 1e4)

	if np.max(redshift) >= z_start:
		raise DarkAgesError("The redshift array is in conflict with the redshift at which to start the evolution of the balck hole. At least one entry in 'redshift' exceeds the value of 'z_start'")

	log10_ini_mass = np.log10(initial_PBH_mass)
	scale = 10**(np.floor(log10_ini_mass)+5)
	initial_PBH_mass *= 1/scale

	# Workaround for dealing with **DarkOptions inside the ODE-solver.
	ODE_to_solve = lambda m,t: PBH_dMdt(m, t, scale=scale, **DarkOptions)
	jac_to_use = lambda m,t: jac(m,t, scale=scale, **DarkOptions)    

	temp_t = 10**np.linspace(np.log10(time_at_z(z_start)), np.log10(time_at_z(1.)), 1e5)
	temp_mass, full_info = solve_ode(ODE_to_solve, initial_PBH_mass, temp_t, Dfun=jac_to_use, full_output=1,mxstep=10000)

	times =  time_at_z(redshift)
	PBH_mass_at_t = np.interp(times, temp_t, temp_mass[:,0])

	out = np.array([redshift, scale*PBH_mass_at_t])
	mask = out[1,:] <= 0.
	out[1,mask] = 0.

	return out
