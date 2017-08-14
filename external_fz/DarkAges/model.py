from .transfer import transfer
from .common import print_info, print_warning, f_function, unscaled, decay_scaling
from .__init__ import DarkAgesError
import numpy as np

class model(object):
	def __init__(self, spec_electrons, spec_photons, normalization, alpha=3):
		self.spec_electrons = spec_electrons
		self.spec_photons = spec_photons
		self.normalization = normalization
		self.alpha_to_use = alpha

	def calc_f(self, transfer_instance):
		if not isinstance(transfer_instance, transfer):
			raise DarkAgesError('You did not include a proper instance of the class "transfer"')
		else:
			red = transfer_instance.z_deposited
			f_func = f_function(transfer_instance.log10E, transfer_instance.z_injected,
                                transfer_instance.z_deposited, self.normalization,
                                transfer_instance.transfer_phot,
                                transfer_instance.transfer_elec,
                                self.spec_photons, self.spec_electrons, alpha=self.alpha_to_use)

			return np.array([red, f_func], dtype=np.float64)

	def save_f(self,transfer_instance, filename):
		f_function = self.calc_f(transfer_instance)
		file_out = open(filename, 'w')
		file_out.write('#z_dep\tf(z)')
		for i in range(len(f_function[0])):
			file_out.write('\n{:.2e}\t{:.4e}'.format(f_function[0,i],f_function[1,i]))
		file_out.close()
		print_info('Saved effective f(z)-curve under "{0}"'.format(filename))

class annihilating_model(model):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m):
		from .__init__ import redshift, logEnergies
		from .common import trapz, unscaled, logConversion
		E = logConversion(logEnergies)
		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec
		#normalization = trapz(tot_spec*E**2, logEnergies)*np.ones_like(redshift)
		normalization = np.ones_like(redshift)*(2*m)
		spec_electrons = np.vectorize(unscaled).__call__(redshift[None,:], ref_el_spec[:,None])
		spec_photons = np.vectorize(unscaled).__call__(redshift[None,:], ref_ph_spec[:,None])
		model.__init__(self, spec_electrons, spec_photons, normalization, 3)

class decaying_model(model):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m, t_dec):
		from .__init__ import redshift, logEnergies
		from .common import trapz, decay_scaling, logConversion
		E = logConversion(logEnergies)
		tot_spec = ref_el_spec + ref_ph_spec + ref_oth_spec
		#normalization = trapz(tot_spec*E**2, logEnergies)*np.ones_like(redshift)
		normalization = np.ones_like(redshift)*(m)
		spec_electrons = np.vectorize(decay_scaling).__call__(redshift[None,:], ref_el_spec[:,None], t_dec)
		spec_photons = np.vectorize(decay_scaling).__call__(redshift[None,:], ref_ph_spec[:,None], t_dec)
		model.__init__(self, spec_electrons, spec_photons, normalization, 0)

class evaporating_model(model):
	def __init__(self, PBH_mass_ini ):
		from .__init__ import redshift, logEnergies
		from .evaporator import PBH_spectrum, PBH_mass_at_z, PBH_dMdt
		from .common import trapz, logConversion, time_at_z
		mass_at_z = PBH_mass_at_z(PBH_mass_ini)
		tmp_mass = np.zeros((mass_at_z.shape[-1]+1), dtype=np.float64)
		tmp_mass[-1] = PBH_mass_ini
		tmp_mass[0:-1] = mass_at_z[-1,:]
		del_mass = np.diff(tmp_mass)
		tmp_z = np.zeros((mass_at_z.shape[-1]+1), dtype=np.float64)
		tmp_z[-1] = 1e4
		tmp_z[0:-1] = mass_at_z[0,:]
		tmp_times = time_at_z(tmp_z)
		del_mass2 = np.vectorize(PBH_dMdt).__call__(mass_at_z[1,:],tmp_times[:-1],1.) * np.diff(tmp_times)

		spec_el = PBH_spectrum( PBH_mass_ini, 'electron')
		spec_ph = PBH_spectrum( PBH_mass_ini, 'gamma')
		spec_all = PBH_spectrum( PBH_mass_ini, 'ALL')
		E = logConversion(logEnergies)
		del_E = np.zeros(redshift.shape, dtype=np.float64)
		for idx in xrange(del_E.shape[0]):
			del_E[idx] = trapz(spec_all[:,idx]*E**2,logEnergies)
			#print('{:.3e} {:.3e} {:.3e} {:.3e}'.format(del_mass[idx], del_mass2[idx], del_E[idx], (del_mass[idx]/del_E[idx])))
		normalization = del_E
		#normalization = del_mass
		#normalization = del_mass2
		#m_ini = PBH_mass_ini * 5.61e23 * 1e9 # mass in eV
		#m_ini = PBH_mass_ini # mass in g
		model.__init__(self, spec_el, spec_ph, normalization, 0)

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
