u"""
.. module:: model
   :synopsis: Definition of the model-class and its derived classes for annihilation, decay, accretion and evaporation of primordial black holes
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

Contains the definition of the base model class :class:`model <DarkAges.model.model>`,
with the basic functions

* :func:`calc_f` to calculate :math:`f(z)`, given an instance of
  :class:`transfer <DarkAges.transfer.transfer>` and
* :func:`model.save_f` to run :func:`calc_f` and saved it in a file.

Also contains derived classes

* :class:`annihilating_model <DarkAges.model.annihilating_model>`
* :class:`annihilating_halos_model <DarkAges.model.annihilating_halos_model>`
* :class:`decaying_model <DarkAges.model.decaying_model>`
* :class:`evaporating_model <DarkAges.model.evaporating_model>`
* :class:`accreting_model <DarkAges.model.accreting_model>`

for the most common energy injection histories.

"""

from __future__ import absolute_import, division, print_function
from builtins import range, object

from .transfer import transfer
from .common import f_function
from .__init__ import DarkAgesError, get_logEnergies, get_redshift, print_info
import numpy as np

class model(object):
	u"""
	Base class to calculate :math:`f(z)` given the injected spectrum
	:math:`\mathrm{d}N / \mathrm{d}E` as a function of *kinetic energy* :math:`E`
	and *redshift* :math:`z+1`
	"""

	def __init__(self, spec_electrons, spec_photons, normalization, logEnergies, alpha=3):
		u"""
		Parameters
		----------
		spec_electrons : :obj:`array-like`
			Array of shape (m,n) containing :math:`\mathrm{d}N / \mathrm{d}E` of
			**electrons** at given redshift :math:`z+1` and
			kinetic energy :math:`E`
		spec_photons : :obj:`array-like`
			Array of shape (m,n) containing :math:`\mathrm{d}N / \mathrm{d}E` of
			**photons** at given redshift :math:`z+1` and
			kinetic energy :math:`E`
		normalization : :obj:`array-like`
			Array of shape (m) with the normalization of the given spectra
			at each given :math:`z_\mathrm{dep.}`.
			(e.g constant array with entries :math:`2m_\mathrm{DM}` for DM-annihilation
			or constant array with entries :math:`m_\mathrm{DM}` for decaying DM)
		alpha : :obj:`int`, :obj:`float`, *optional*
			Exponent to specify the comoving scaling of the
			injected spectra.
			(3 for annihilation and 0 for decaying species
			`c.f. ArXiv1801.01871 <https://arxiv.org/abs/1801.01871>`_).
			If not specified annihilation is assumed.
		"""

		self.logEnergies = logEnergies
		self.spec_electrons = spec_electrons
		self.spec_photons = spec_photons
		self.normalization = normalization
		self.alpha_to_use = alpha

	def calc_f(self, transfer_instance, **DarkOptions):
		u"""Returns :math:`f(z)` for a given set of transfer functions
		:math:`T(z_{dep}, E, z_{inj})`

		Parameters
		----------
		transfer_instance : :obj:`class`
			Initialized instace of :class:`transfer <DarkAges.transfer.transfer>`

		Returns
		-------
		:obj:`array-like`
			Array (:code:`shape=(2,n)`) containing :math:`z_\mathrm{dep}+1` in the first column
			and :math:`f(z_\mathrm{dep})` in the second column.
		"""

		if not isinstance(transfer_instance, transfer):
			raise DarkAgesError('You did not include a proper instance of the class "transfer"')
		else:
			red = transfer_instance.z_deposited

			f_func = f_function(transfer_instance.log10E,self.logEnergies, transfer_instance.z_injected,
                                transfer_instance.z_deposited, self.normalization,
                                transfer_instance.transfer_phot,
                                transfer_instance.transfer_elec,
                                self.spec_photons, self.spec_electrons, alpha=self.alpha_to_use, **DarkOptions)

			return np.array([red, f_func], dtype=np.float64)

	def save_f(self,transfer_instance, filename, **DarkOptions):
		u"""Saves the table :math:`z_\mathrm{dep.}`, :math:`f(z_\mathrm{dep})` for
		a given set of transfer functions :math:`T(z_{dep}, E, z_{inj})` in a file.

		Parameters
		----------
		transfer_instance : :obj:`class`
			Initialized instace of :class:`transfer <DarkAges.transfer.transfer>`
		filename : :obj:`str`
			Self-explanatory
		"""

		f_function = self.calc_f(transfer_instance,**DarkOptions)
		file_out = open(filename, 'w')
		file_out.write('#z_dep\tf(z)')
		for i in range(len(f_function[0])):
			file_out.write('\n{:.2e}\t{:.4e}'.format(f_function[0,i],f_function[1,i]))
		file_out.close()
		print_info('Saved effective f(z)-curve under "{0}"'.format(filename))

class annihilating_model(model):
	u"""Derived instance of the class :class:`model <DarkAges.model.model>` for the case of an annihilating
	species.

	Inherits all methods of :class:`model <DarkAges.model.model>`
	"""

	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,logEnergies = None,redshift=None, **DarkOptions):
		u"""
		At initialization the reference spectra are read and the double-differential
		spectrum :math:`\\frac{\\mathrm{d}^2 N(t,E)}{\\mathrm{d}E\\mathrm{d}t}` needed for
		the initialization inherited from :class:`model <DarkAges.model.model>` is calculated by

		.. math::
			\\frac{\\mathrm{d}^2 N(t,E)}{\\mathrm{d}E\\mathrm{d}t} = C \\cdot\\frac{\\mathrm{d}N(E)}{\\mathrm{d}E}

		where :math:`C` is a constant independent of :math:`t` (:math:`z`) and :math:`E`

		Parameters
		----------
		ref_el_spec : :obj:`array-like`
			Reference spectrum (:code:`shape = (k,l)`) :math:`\mathrm{d}N / \mathrm{d}E` of **electrons**
		ref_ph_spec : :obj:`array-like`
			Reference spectrum (:code:`shape = (k,l)`) :math:`\mathrm{d}N / \mathrm{d}E` of **photons**
		ref_oth_spec : :obj:`array-like`
			Reference spectrum (:code:`shape = (k,l)`) :math:`\mathrm{d}N / \mathrm{d}E` of particles
			not interacting with the erly IGM (e.g. **protons** and **neutrinos**).
			This is neede for the proper normalization of the electron- and photon-spectra.
		m : :obj:`float`
			Mass of the DM-candidate (*in units of* :math:`\\mathrm{GeV}`)
		logEnergies : :obj:`array-like`, optional
			Array (:code:`shape = (l)`) of the logarithms of the kinetic energies of the particles
			(*in units of* :math:`\\mathrm{eV}`) to the base 10.
			If not specified, the standard array provided by
			:class:`the initializer <DarkAges.__init__>`  is taken.
		redshift : :obj:`array-like`, optional
			Array (:code:`shape = (k)`) with the values of :math:`z+1`. Used for
			the calculation of the double-differential spectra.
			If not specified, the standard array provided by
			:mod:`the initializer <DarkAges.__init__>`  is taken.
		"""

		if logEnergies is None:
			logEnergies = get_logEnergies()
		if redshift is None:
			redshift = get_redshift()

		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec

		norm_by = DarkOptions.get('normalize_spectrum_by','energy_integral')
		if norm_by == 'energy_integral':
			from .common import trapz, logConversion
			E = logConversion(logEnergies)
			if len(E) > 1:
				normalization = trapz(tot_spec*E**2*np.log(10), logEnergies)*np.ones_like(redshift)
			else:
				normalization = (tot_spec*E)[0]
		elif norm_by == 'mass':
			normalization = np.ones_like(redshift)*(2*m)
		else:
			raise DarkAgesError('I did not understand your input of "normalize_spectrum_by" ( = {:s}). Please choose either "mass" or "energy_integral"'.format(norm_by))

		spec_electrons = np.zeros((len(tot_spec),len(redshift)))
		spec_photons = np.zeros((len(tot_spec),len(redshift)))
		spec_electrons[:,:] = ref_el_spec[:,None]
		spec_photons[:,:] = ref_ph_spec[:,None]

		model.__init__(self, spec_electrons, spec_photons, normalization,logEnergies, 3)

class annihilating_halos_model(model):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,zh,fh,logEnergies=None,redshift=None, **DarkOptions):

		from .special_functions import boost_factor_halos

		def scaling_boost_factor(redshift,spec_point,zh,fh):
			ret = spec_point*boost_factor_halos(redshift,zh,fh)
			return ret

		if logEnergies is None:
			logEnergies = get_logEnergies()
		if redshift is None:
			redshift = get_redshift()

		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec

		norm_by = DarkOptions.get('normalize_spectrum_by','energy_integral')
		if norm_by == 'energy_integral':
			from .common import trapz, logConversion
			E = logConversion(logEnergies)
			if len(E) > 1:
				normalization = trapz(tot_spec*E**2*np.log(10), logEnergies)*np.ones_like(redshift)
			else:
				normalization = (tot_spec*E)[0]
		elif norm_by == 'mass':
			normalization = np.ones_like(redshift)*(2*m)
		else:
			raise DarkAgesError('I did not understand your input of "normalize_spectrum_by" ( = {:s}). Please choose either "mass" or "energy_integral"'.format(norm_by))
		normalization /= boost_factor_halos(redshift,zh,fh)

		spec_electrons = np.vectorize(scaling_boost_factor).__call__(redshift[None,:],ref_el_spec[:,None],zh,fh)
		spec_photons = np.vectorize(scaling_boost_factor).__call__(redshift[None,:],ref_ph_spec[:,None],zh,fh)

		model.__init__(self, spec_electrons, spec_photons, normalization, logEnergies,3)


class decaying_model(model):
	u"""Derived instance of the class :class:`model <DarkAges.model.model>` for the case of a decaying
	species.

	Inherits all methods of :class:`model <DarkAges.model.model>`
	"""

	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,t_dec,logEnergies=None,redshift=None, **DarkOptions):
		u"""At initialization the reference spectra are read and the double-differential
		spectrum :math:`\\frac{\\mathrm{d}^2 N(t,E)}{\\mathrm{d}E\\mathrm{d}t}` needed for
		the initialization inherited from :class:`model <DarkAges.model.model>` is calculated by

		.. math::
			\\frac{\\mathrm{d}^2 N(t,E)}{\\mathrm{d}E\\mathrm{d}t} = C \\cdot\\exp{\\left(\\frac{-t(z)}{\\tau}\\right)} \\cdot \\frac{\\mathrm{d}N(E)}{\\mathrm{d}E}

		where :math:`C` is a constant independent of :math:`t` (:math:`z`) and :math:`E`

		Parameters
		----------
		ref_el_spec : :obj:`array-like`
			Reference spectrum (:code:`shape = (k,l)`) :math:`\mathrm{d}N / \mathrm{d}E` of **electrons**
		ref_ph_spec : :obj:`array-like`
			Reference spectrum (:code:`shape = (k,l)`) :math:`\mathrm{d}N / \mathrm{d}E` of **photons**
		ref_oth_spec : :obj:`array-like`
			Reference spectrum (:code:`shape = (k,l)`) :math:`\mathrm{d}N / \mathrm{d}E` of particles
			not interacting with the early IGM (e.g. **protons** and **neutrinos**).
			This is needed for the proper normalization of the electron- and photon-spectra.
		m : :obj:`float`
			Mass of the DM-candidate (*in units of* :math:`\\mathrm{GeV}`)
		t_dec : :obj:`float`
			Lifetime (Time after which the number of particles dropped down to
			a factor of :math:`1/e`) of the DM-candidate
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

		def _decay_scaling(redshift, spec_point, lifetime):
			from .common import time_at_z
			ret = spec_point*np.exp(-time_at_z(redshift) / lifetime)
			return ret

		if logEnergies is None:
			logEnergies = get_logEnergies()
		if redshift is None:
			redshift = get_redshift()

		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec

		norm_by = DarkOptions.get('normalize_spectrum_by','energy_integral')
		if norm_by == 'energy_integral':
			from .common import trapz, logConversion
			E = logConversion(logEnergies)
			if len(E) > 1:
				normalization = trapz(tot_spec*E**2*np.log(10), logEnergies)*np.ones_like(redshift)
			else:
				normalization = (tot_spec*E)[0]
		elif norm_by == 'mass':
			normalization = np.ones_like(redshift)*(m)
		else:
			raise DarkAgesError('I did not understand your input of "normalize_spectrum_by" ( = {:s}). Please choose either "mass" or "energy_integral"'.format(norm_by))

		spec_electrons = np.vectorize(_decay_scaling).__call__(redshift[None,:], ref_el_spec[:,None], t_dec)
		spec_photons = np.vectorize(_decay_scaling).__call__(redshift[None,:], ref_ph_spec[:,None], t_dec)

		model.__init__(self, spec_electrons, spec_photons, normalization, logEnergies,0)

class evaporating_model(model):
	u"""Derived instance of the class :class:`model <DarkAges.model.model>` for the case of evaporating
	primordial black holes (PBH) as a candidate of DM

	Inherits all methods of :class:`model <DarkAges.model.model>`
	"""

	def __init__(self, PBH_mass_ini, logEnergies=None, redshift=None, **DarkOptions):
		u"""
		At initialization evolution of the PBH mass is calculated with
		:func:`PBH_mass_at_z <DarkAges.evaporator.PBH_mass_at_z>` and the
		double-differential spectrum :math:`\mathrm{d}^2 N(z,E) / \mathrm{d}E\mathrm{d}z`
		needed for the initialization inherited from :class:`model <DarkAges.model.model>` is calculated
		according to :func:`PBH_spectrum <DarkAges.evaporator.PBH_spectrum>`

		Parameters
		----------
		PBH_mass_ini : :obj:`float`
			Initial mass of the primordial black hole (*in units of* :math:`\\mathrm{g}`)
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

		from .evaporator import PBH_spectrum_at_m, PBH_mass_at_z, PBH_dMdt
		from .common import trapz, logConversion, time_at_z, nan_clean

		include_secondaries=DarkOptions.get('PBH_with_secondaries',True)

		if logEnergies is None:
			logEnergies = get_logEnergies()
		if redshift is None:
			redshift = get_redshift()

		mass_at_z = PBH_mass_at_z(PBH_mass_ini, redshift=redshift, **DarkOptions)
		dMdt_at_z = (-1)*np.vectorize(PBH_dMdt).__call__(mass_at_z[-1,:],np.ones_like(mass_at_z[0,:]))

		E = logConversion(logEnergies)
		E_sec = 1e-9*E
		E_prim = 1e-9*E

		# Total spectrum (for normalization)
		spec_all = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'ALL', **DarkOptions)
		del_E = np.zeros(redshift.shape, dtype=np.float64)
		corr = np.ones(redshift.shape, dtype=np.float64)
		for idx in range(del_E.shape[0]):
			del_E[idx] = trapz(spec_all[:,idx]*E**2*np.log(10),(logEnergies))
			if del_E[idx] > 0.0:
				corr[idx] = dMdt_at_z[idx]/del_E[idx]
		normalization = del_E*corr

		# Primary spectra
		prim_spec_el = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'electron', **DarkOptions)*corr[None,:]
		prim_spec_ph = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'gamma', **DarkOptions)*corr[None,:]
		prim_spec_muon = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'muon', **DarkOptions)*corr[None,:]
		prim_spec_pi0 = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'pi0', **DarkOptions)*corr[None,:]
		prim_spec_piCh = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'piCh', **DarkOptions)*corr[None,:]

		# full spectra (including secondaries)
		if include_secondaries:
			from .special_functions import secondaries_from_simple_decay
			sec_from_pi0 = secondaries_from_simple_decay(E_sec[:,None],E_prim[None,:],'pi0')
			sec_from_pi0 /= (trapz(np.sum(sec_from_pi0, axis=2),E_sec,axis=0))[None,:,None]
			sec_from_piCh = secondaries_from_simple_decay(E_sec[:,None],E_prim[None,:],'piCh')
			sec_from_piCh /= (trapz(np.sum(sec_from_piCh, axis=2),E_sec,axis=0))[None,:,None]
			sec_from_muon = secondaries_from_simple_decay(E_sec[:,None],E_prim[None,:],'muon')
			sec_from_muon /= (trapz(np.sum(sec_from_muon, axis=2),E_sec,axis=0))[None,:,None]
		else:
			sec_from_pi0 = np.zeros((len(E_sec),len(E_prim),3), dtype=np.float64)
			sec_from_piCh = np.zeros((len(E_sec),len(E_prim),3), dtype=np.float64)
			sec_from_muon = np.zeros((len(E_sec),len(E_prim),3), dtype=np.float64)

		spec_el = np.zeros_like(prim_spec_el)
		spec_el += prim_spec_el
		spec_el += trapz((sec_from_pi0[:,:,None,0])*prim_spec_pi0[None,:,:],E_prim,axis=1)
		spec_el += trapz((sec_from_piCh[:,:,None,0])*prim_spec_piCh[None,:,:],E_prim,axis=1)
		spec_el += trapz((sec_from_muon[:,:,None,0])*prim_spec_muon[None,:,:],E_prim,axis=1)
		spec_el =  nan_clean(spec_el)

		spec_ph = np.zeros_like(prim_spec_ph)
		spec_ph += prim_spec_ph
		spec_ph += trapz((sec_from_pi0[:,:,None,1])*prim_spec_pi0[None,:,:],E_prim,axis=1)
		spec_ph += trapz((sec_from_piCh[:,:,None,1])*prim_spec_piCh[None,:,:],E_prim,axis=1)
		spec_ph += trapz((sec_from_muon[:,:,None,1])*prim_spec_muon[None,:,:],E_prim,axis=1)
		spec_ph = nan_clean(spec_ph)

		model.__init__(self, spec_el, spec_ph, normalization, logEnergies,0)

class accreting_model(model):
	u"""Derived instance of the class :class:`model <DarkAges.model.model>` for
	the case of accreting primordial black holes (PBH) as a candidate of DM.

	Inherits all methods of :class:`model <DarkAges.model.model>`
	"""

	def __init__(self, PBH_mass, recipe, logEnergies=None, redshift=None, **DarkOptions):
		u"""At initialization the reference spectra are read and the luminosity
		spectrum :math:`L_{\\omega}` needed for the initialization inherited
		from :class:`model <DarkAges.model.model>` is calculated by

		.. math::
			L_{\\omega} = \\Theta(\\omega -\\omega_\\mathrm{min})w^{-a}\\exp(-\\frac{\\omega}{T_s})

		where :math:`T_s\\simeq 200\\,\\mathrm{keV}`, :math:`a=-2.5+\\frac{\\log(M)}{3}` and
		:math:`\\omega_\\mathrm{min} = \\left(\\frac{10}{M}\\right)^{\\frac{1}{2}}`
		if :code:`recipe = disk_accretion` or

		..  math::
			L_\\omega = w^{-a}\\exp(-\\frac{\\omega}{T_s})

		where :math:`T_s\\simeq 200\\,\\mathrm{keV}` if
		:code:`recipe = spherical_accretion`.

		Parameters
		----------
		PBH_mass : :obj:`float`
			Mass of the primordial black hole (*in units of* :math:`M_\\odot`)
		recipe : :obj:`string`
			Recipe setting the luminosity and the rate of the accretion
			(`spherical_accretion` taken from 1612.05644 and `disk_accretion`
			from 1707.04206)
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

		if logEnergies is None:
			logEnergies = get_logEnergies()
		if redshift is None:
			redshift = get_redshift()

		from .common import trapz,  logConversion
		from .special_functions import luminosity_accreting_bh
		E = logConversion(logEnergies)
		spec_ph = luminosity_accreting_bh(E,recipe,PBH_mass)
		spec_el = np.zeros_like(spec_ph)
		spec_oth = np.zeros_like(spec_ph)
		normalization = trapz((spec_ph+spec_el)*E**2*np.log(10),logEnergies)*np.ones_like(redshift)

		spec_photons = np.zeros((len(spec_el),len(redshift)))
		spec_photons[:,:] = spec_ph[:,None]
		spec_electrons = np.zeros((len(spec_el),len(redshift)))

		model.__init__(self, spec_electrons, spec_photons, normalization, logEnergies, 0)
