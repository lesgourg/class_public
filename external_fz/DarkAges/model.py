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
* :class:`decaying_model <DarkAges.model.decaying_model>`
* :class:`evaporating_model <DarkAges.model.evaporating_model>`
* :class:`accreting_model <DarkAges.model.accreting_model>`

for the most common energy injection histories.

"""
from .transfer import transfer
from .common import print_info, print_warning, f_function, boost_factor_halos, scaling_boost_factor
from .__init__ import DarkAgesError
import numpy as np

class model(object):
	u"""
	Base class to calculate :math:`f(z)` given the injected spectrum
	:math:`\mathrm{d}N / \mathrm{d}E` as a function of *kinetic energy* :math:`E`
	and *redshift* :math:`z+1`
	"""

	def __init__(self, spec_electrons, spec_photons, normalization, alpha=3):
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
			`c.f. ArXivXXXX.YYYY <https://arxiv.org/abs/XXXX.YYYY>`_).
			If not specified annihilation is assumed.
		"""

		self.spec_electrons = spec_electrons
		self.spec_photons = spec_photons
		self.normalization = normalization
		self.alpha_to_use = alpha

	def calc_f(self, transfer_instance):
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
            # interpolate_transfer(transfer_instance.transfer_phot,transfer_instance.transfer_elec,transfer_instance.log10E, transfer_instance.z_injected,
            #                     transfer_instance.z_deposited)
			red = transfer_instance.z_deposited
			f_func = f_function(transfer_instance.log10E, transfer_instance.z_injected,
                                transfer_instance.z_deposited, self.normalization,
                                transfer_instance.transfer_phot,
                                transfer_instance.transfer_elec,
                                self.spec_photons, self.spec_electrons, alpha=self.alpha_to_use)

			return np.array([red, f_func], dtype=np.float64)

	def save_f(self,transfer_instance, filename):
		u"""Saves the table :math:`z_\mathrm{dep.}`, :math:`f(z_\mathrm{dep})` for
		a given set of transfer functions :math:`T(z_{dep}, E, z_{inj})` in a file.

		Parameters
		----------
		transfer_instance : :obj:`class`
			Initialized instace of :class:`transfer <DarkAges.transfer.transfer>`
		filename : :obj:`str`
			Self-explanatory
		"""

		f_function = self.calc_f(transfer_instance)
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

	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,logEnergies=None,redshift=None):
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

		def _unscaled(redshift, spec_point):
			ret = spec_point*np.ones_like(redshift)
			return ret

		if logEnergies is None:
			from .__init__ import logEnergies as cheese
			logEnergies = cheese
		if redshift is None:
			from .__init__ import redshift as ham
			redshift = ham
		from .common import trapz, logConversion

		E = logConversion(logEnergies)
		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec
		#normalization = trapz(tot_spec*E**2, logEnergies)*np.ones_like(redshift)
		normalization = np.ones_like(redshift)*(2*m)
		spec_electrons = np.vectorize(_unscaled).__call__(redshift[None,:], ref_el_spec[:,None])
		spec_photons = np.vectorize(_unscaled).__call__(redshift[None,:], ref_ph_spec[:,None])
		model.__init__(self, spec_electrons, spec_photons, normalization, 3)

class annihilating_halos_model(model):
        def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,zh,fh):
                from .__init__ import redshift, logEnergies
                from .common import trapz, unscaled, logConversion
                E = logConversion(logEnergies)
                tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec
                #normalization = trapz(tot_spec*E**2, logEnergies)*np.ones_like(redshift)
                normalization = np.ones_like(redshift)*(2*m)/boost_factor_halos(redshift,zh,fh)
                spec_electrons = np.vectorize(scaling_boost_factor).__call__(redshift[None,:],ref_el_spec[:,None],zh,fh)
                spec_photons = np.vectorize(scaling_boost_factor).__call__(redshift[None,:],ref_ph_spec[:,None],zh,fh)
                model.__init__(self, spec_electrons, spec_photons, normalization, 3)


class decaying_model(model):
	u"""Derived instance of the class :class:`model <DarkAges.model.model>` for the case of a decaying
	species.

	Inherits all methods of :class:`model <DarkAges.model.model>`
	"""

	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,t_dec,logEnergies=None,redshift=None):
		u"""
		At initialization the reference spectra are read and the double-differential
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
			from .__init__ import logEnergies as cheese
			logEnergies = cheese
		if redshift is None:
			from .__init__ import redshift as ham
			redshift = ham
		from .common import trapz, logConversion

		E = logConversion(logEnergies)
		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec
		#normalization = trapz(tot_spec*E**2, logEnergies)*np.ones_like(redshift)
		normalization = np.ones_like(redshift)*(m)
		spec_electrons = np.vectorize(_decay_scaling).__call__(redshift[None,:], ref_el_spec[:,None], t_dec)
		spec_photons = np.vectorize(_decay_scaling).__call__(redshift[None,:], ref_ph_spec[:,None], t_dec)
		model.__init__(self, spec_electrons, spec_photons, normalization, 0)

class evaporating_model(model):
    u"""Derived instance of the class :class:`model <DarkAges.model.model>` for the case of evaporating
    primordial black holes (PBH) as a candidate of DM

    Inherits all methods of :class:`model <DarkAges.model.model>`
    """

    def __init__(self, PBH_mass_ini, logEnergies=None, redshift=None):
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

        if logEnergies is None:
            from .__init__ import logEnergies as cheese
            logEnergies = cheese
        if redshift is None:
            from .__init__ import redshift as ham
            redshift = ham

        mass_at_z = PBH_mass_at_z(PBH_mass_ini, redshift=redshift)
        E = logConversion(logEnergies)
        E_sec = 1e-9*E
        E_prim = 1e-9*E

        # Primary spectra
        prim_spec_el = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'electron')
        prim_spec_ph = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'gamma')
        prim_spec_muon = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'muon')
        prim_spec_pi0 = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'pi0')
        prim_spec_piCh = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'piCh')

        # full spectra (including secondaries)
        from .special_functions import secondaries_from_pi0, secondaries_from_piCh, secondaries_from_muon
        sec_from_pi0 = secondaries_from_pi0(E_sec[:,None],E_prim[None,:])
        #sec_from_pi0 /= (trapz(np.sum(sec_from_pi0, axis=2)*E_sec[:,None],E_sec,axis=0)/(E_prim))[None,:,None]
        #sec_from_pi0 /= (trapz(np.sum(sec_from_pi0, axis=2),E_sec,axis=0))[None,:,None]
        sec_from_piCh = secondaries_from_piCh(E_sec[:,None],E_prim[None,:])
        #sec_from_piCh /= (trapz(np.sum(sec_from_piCh, axis=2)*E_sec[:,None],E_sec,axis=0)/(E_prim))[None,:,None]
        #sec_from_piCh /= (trapz(np.sum(sec_from_piCh, axis=2),E_sec,axis=0))[None,:,None]
        sec_from_muon = secondaries_from_muon(E_sec[:,None],E_prim[None,:])
        #sec_from_muon /= (trapz(np.sum(sec_from_muon, axis=2)*E_sec[:,None],E_sec,axis=0)/(E_prim))[None,:,None]
        #sec_from_muon /= (trapz(np.sum(sec_from_muon, axis=2),E_sec,axis=0))[None,:,None]

        #convol_norm_pi0 = trapz(np.sum(sec_from_pi0, axis=2)*E_prim[None,:],E_prim,axis=1)
        #convol_norm_piCh = trapz(np.sum(sec_from_piCh, axis=2)*E_prim[None,:],E_prim,axis=1)
        #convol_norm_muon = trapz(np.sum(sec_from_muon, axis=2)*E_prim[None,:],E_prim,axis=1)

        convol_norm_pi0 = np.ones((1,))
        convol_norm_piCh = np.ones((1,))
        convol_norm_muon = np.ones((1,))

        spec_el = np.zeros_like(prim_spec_el)
        spec_el += prim_spec_el
        spec_el += trapz((sec_from_pi0[:,:,None,0])*prim_spec_pi0[None,:,:],E_prim,axis=1)/convol_norm_pi0[:,None]
        spec_el += trapz((sec_from_piCh[:,:,None,0])*prim_spec_piCh[None,:,:],E_prim,axis=1)/convol_norm_piCh[:,None]
        spec_el += trapz((sec_from_muon[:,:,None,0])*prim_spec_muon[None,:,:],E_prim,axis=1)/convol_norm_muon[:,None]
        spec_el =  nan_clean(spec_el)

        spec_ph = np.zeros_like(prim_spec_ph)
        spec_ph += prim_spec_ph
        spec_ph += trapz((sec_from_pi0[:,:,None,1])*prim_spec_pi0[None,:,:],E_prim,axis=1)/convol_norm_pi0[:,None]
        spec_ph += trapz((sec_from_piCh[:,:,None,1])*prim_spec_piCh[None,:,:],E_prim,axis=1)/convol_norm_piCh[:,None]
        spec_ph += trapz((sec_from_muon[:,:,None,1])*prim_spec_muon[None,:,:],E_prim,axis=1)/convol_norm_muon[:,None]
        spec_ph = nan_clean(spec_ph)

        # Total spectrum (for normalization)
        spec_all = PBH_spectrum_at_m( mass_at_z[-1,:], logEnergies, 'ALL')
        del_E = np.zeros(redshift.shape, dtype=np.float64)
        for idx in xrange(del_E.shape[0]):
            del_E[idx] = trapz(spec_all[:,idx]*E**2,(logEnergies))
            normalization = del_E

        model.__init__(self, spec_el, spec_ph, normalization, 0)

class accreting_model(model):
    u"""Derived instance of the class :class:`model <DarkAges.model.model>` for the case of accreting
    primordial black holes (PBH) as a candidate of DM

    Inherits all methods of :class:`model <DarkAges.model.model>`
    """

    def __init__(self, PBH_mass, recipe, logEnergies=None, redshift=None):
        u"""
        At initialization the reference spectra are read and the luminosity
        spectrum :math:`L_\omega` needed for
        the initialization inherited from :class:`model <DarkAges.model.model>` is calculated by

        .. math::
        	L_\omega = \Theta(\omega - \omega_{\rm min})w^{-a}\exp(-\omega/T_s)
        where :math:`T_s\simeq 200 keV`, :math:`a=-2.5+\log(M)/3` and :math:`\omega_{\rm min} = (10/M)^{1/2}` if recipe = disk_accretion or
        .. math::
            L_\omega = w^{-a}\exp(-\omega/T_s)
        where :math:`T_s\simeq 200 keV` if recipe = spherical_accretion.

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

        from .__init__ import redshift, logEnergies
        from .common import trapz,  logConversion
        from .special_functions import luminosity_accreting_bh
        E = logConversion(logEnergies)
        spec_ph = luminosity_accreting_bh(E,recipe,PBH_mass)
        spec_el = np.zeros_like(spec_ph)
        def _unscaled(redshift, spec_point):
			ret = spec_point*np.ones_like(redshift)
			return ret

        spec_electrons = np.vectorize(_unscaled).__call__(redshift[None,:], spec_el[:,None])
        spec_photons = np.vectorize(_unscaled).__call__(redshift[None,:], spec_ph[:,None])
        normalization = trapz((spec_ph+spec_el)*E,logEnergies)*np.ones_like(redshift)
        # print normalization, spec_photons
        model.__init__(self, spec_electrons, spec_photons, normalization, 0)



#### OLD ######

'''
class old_model(object):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,history='annihilation', t_dec = None):
		self._is_initialized = False
		self.mass = m
		self.injection_hist = history
		self.photon_spec = ref_ph_spec
		self.electron_spec = ref_el_spec
		self.total_spec = ref_ph_spec + ref_el_spec + ref_oth_spec
		if history == 'decay':
			self.decay_time = t_dec

	def make_z_dependent_spectrum(self, redshift, ref_spectrum, z_dependence=unscaled, *args_for_z_dependence):
		out_spec = np.zeros(shape=(len(ref_spectrum),len(redshift)), dtype=np.float64)
		for idx_E in xrange(out_spec.shape[0]):
			spec_point = ref_spectrum[idx_E]
			out_spec[idx_E,:] = z_dependence(redshift, spec_point, *args_for_z_dependence)
		return out_spec

	def get_normalization(self, redshift):
		out_norm = np.ones_like(redshift)
		if self.injection_hist == 'annihilation' or self.injection_hist == 'PBH':
			out_norm *= 2*self.mass
		elif self.injection_hist == 'decay':
			out_norm *= self.mass
		return out_norm

	def calc_f(self, transfer_instance):
		alpha_dict = {'annihilation':3, 'decay':0, 'PBH':0}
		if not isinstance(transfer_instance, transfer):
			print_warning('You did not include a proper instance of the class "transfer"')
			return -1
		else:
			red = transfer_instance.z_deposited
			if not self._is_initialized:
				self._is_initialized = True
				if self.injection_hist == 'annihilation' or self.injection_hist == 'PBH':
					self.z_dependent_electrons = self.make_z_dependent_spectrum(red, self.electron_spec, unscaled)
					self.z_dependent_photons = self.make_z_dependent_spectrum(red, self.photon_spec, unscaled)
				elif self.injection_hist == 'decay':
					self.z_dependent_electrons = self.make_z_dependent_spectrum(red, self.electron_spec, decay_scaling, self.decay_time)
					self.z_dependent_photons = self.make_z_dependent_spectrum(red, self.photon_spec, decay_scaling, self.decay_time)
			if self.injection_hist in alpha_dict:
				alpha_to_use = alpha_dict[self.injection_hist]
				self.normalization = self.get_normalization(red)
			else:
				raise DarkAgesError('The code can not deal with the injection history >> {0} << (yet)'.format(self.injection_hist))

			f_func = f_function(transfer_instance.log10E, transfer_instance.z_injected,
                                transfer_instance.z_deposited, self.normalization,
                                transfer_instance.transfer_phot,
                                transfer_instance.transfer_elec,
                                self.z_dependent_photons, self.z_dependent_electrons, alpha=alpha_to_use)

			return np.array([red, f_func], dtype=np.float64)

	def save_f(self,transfer_instance, filename):
		f_function = self.calc_f(transfer_instance)
		file_out = open(filename, 'w')
		file_out.write('#z_dep\tf(z)')
		for i in range(len(f_function[0])):
			file_out.write('\n{:.2e}\t{:.4e}'.format(f_function[0,i],f_function[1,i]))
		file_out.close()
		print_info('Saved effective f(z)-curve under "{0}"'.format(filename))

class annihilating_model2(old_model):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,):
		old_model.__init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,history='annihilation')

class decaying_model2(old_model):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m, t_dec):
		old_model.__init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,history='decay', t_dec=t_dec)
'''
