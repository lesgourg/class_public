import numpy as np
from .common import print_warning

class transfer(object):
	def __init__(self, infile):
		#print 'Initializing the transfer functions'
		data = np.genfromtxt(infile, unpack=True, usecols=(0,1,2,3,4), dtype=np.float64 )
		self.z_injected = np.unique(data[2]).astype(np.float64)
		self.z_deposited = np.unique(data[0]).astype(np.float64)
		self.log10E = np.unique(data[1]).astype(np.float64)
		l1 = len(self.z_deposited)
		l2 = len(self.log10E)
		l3 = len(self.z_injected)
		self.transfer_phot = data[4].reshape(l1,l2,l3).astype(np.float64)
		self.transfer_elec = data[3].reshape(l1,l2,l3).astype(np.float64)

def transfer_dump(transfer_instance, outfile):
	import dill
	if not isinstance(transfer_instance, transfer):
		print_warning('You did not include a proper instance of the class "transfer"')
		return -1
	with open(outfile, 'wb') as f_dump:
		dill.dump(transfer_instance, f_dump)
	return 0

def transfer_load(infile):
	import dill
	loaded_transfer = dill.load(open(infile, 'rb'))
	if not isinstance(loaded_transfer, transfer):
		print_warning('The file {} does not provide a proper instance of the class "transfer"'.format(infile))
		return -1
	else:
		return loaded_transfer
