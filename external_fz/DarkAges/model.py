from .transfer import transfer
from .common import print_info, print_warning, print_error, f_function, unscaled, decay_scaling
import numpy as np

class model(object):
	def __init__(self,ref_el_spec,ref_ph_spec,ref_oth_spec,m,history='annihilation', t_dec = None):
		self.z_dependence_initialized = False
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

	def calc_f(self, transfer_instance):
		alpha_dict = {'annihilation':3, 'decay':0}
		if not isinstance(transfer_instance, transfer):
			print_warning('You did not include a proper instance of the class "transfer"')
			return -1
		else:
			red = transfer_instance.z_deposited
			if not self.z_dependence_initialized:
				self.z_dependence_initialized = True
				if self.injection_hist == 'annihilation':
					self.z_dependent_electrons = self.make_z_dependent_spectrum(red, self.electron_spec, unscaled)
					self.z_dependent_photons = self.make_z_dependent_spectrum(red, self.photon_spec, unscaled)
				elif self.injection_hist == 'decay':
					self.z_dependent_electrons = self.make_z_dependent_spectrum(red, self.electron_spec, decay_scaling, self.decay_time)
					self.z_dependent_photons = self.make_z_dependent_spectrum(red, self.photon_spec, decay_scaling, self.decay_time)
			if self.injection_hist in alpha_dict:
				alpha_to_use = alpha_dict[self.injection_hist]
			else:
				print_error('The code can not deal with the injection history >> {} << (yet)'.format(self.injection_hist))
			f_func = f_function(transfer_instance.log10E, transfer_instance.z_injected,
                                transfer_instance.z_deposited, self.mass,
                                transfer_instance.transfer_phot,
                                transfer_instance.transfer_elec,
                                self.z_dependent_photons, self.z_dependent_electrons, self.total_spec, alpha=alpha_to_use)

			return np.array([red, f_func], dtype=np.float64)

	def save_f(self,transfer_instance, filename):
		f_function = self.calc_f(transfer_instance)
		file_out = open(filename, 'w')
		file_out.write('#z_dep\tf(z)')
		for i in range(len(f_function[0])):
			file_out.write('\n{:.2e}\t{:.4e}'.format(f_function[0,i],f_function[1,i]))
		file_out.close()
		print_info('Saved effective f(z)-curve under "{}"'.format(filename))
