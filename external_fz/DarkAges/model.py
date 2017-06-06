from .transfer import transfer
from .common import print_info, print_warning, print_error, f_function_annihilation, f_function_decay
import numpy as np

class model(object):
	def __init__(self,el_spec,ph_spec,oth_spec,m,history='annihilation', t_dec = None):
		self.mass = m
		self.photon_spec = ph_spec
		self.electron_spec = el_spec
		self.total_spec = ph_spec+el_spec+oth_spec
		self.injection_hist = history
		if history == 'decay':
			self.decay_time = t_dec

	def calc_f(self, transfer_instance):
		if not isinstance(transfer_instance, transfer):
			print_warning('You did not include a proper instance of the class "transfer"')
			return -1
		else:
			red = transfer_instance.z_deposited
			if self.injection_hist == 'annihilation':
				f_func = f_function_annihilation(transfer_instance.log10E, transfer_instance.z_injected,
                                                 transfer_instance.z_deposited, self.mass,                                                  transfer_instance.transfer_phot,
                                                 transfer_instance.transfer_elec,
                                                 self.photon_spec, self.electron_spec, self.total_spec)
			elif self.injection_hist == 'decay':
				f_func = f_function_decay(transfer_instance.log10E, transfer_instance.z_injected,
                                          transfer_instance.z_deposited, self.mass, self.decay_time,
                                          transfer_instance.transfer_phot,
                                          transfer_instance.transfer_elec,
                                          self.photon_spec, self.electron_spec, self.total_spec)
			else:
				print_error('The code can not deal with the injection history >> {} << (yet)'.format(self.injection_hist))
			return np.array([red, f_func], dtype=np.float64)

	def save_f(self,transfer_instance, filename):
		f_function = self.calc_f(transfer_instance)
		file_out = open(filename, 'w')
		file_out.write('#z_dep\tf(z)')
		for i in range(len(f_function[0])):
			file_out.write('\n{:.2e}\t{:.4e}'.format(f_function[0,i],f_function[1,i]))
		file_out.close()
		print_info('Saved effective f(z)-curve under "{}"'.format(filename))
