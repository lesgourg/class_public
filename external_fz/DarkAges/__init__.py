import numpy as np
import sys
import os
os.environ['DARKAGES_BASE'] = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.insert(0, os.environ['DARKAGES_BASE'])

if os.getppid() == int(os.environ['DARKAGES_TOPLEVEL_PID']):
	import dill
	with open(os.path.join( os.environ['DARKAGES_BASE'], 'pid_{:d}'.format(os.getppid()) ), 'rb') as options_dumper:
		options = dill.load(options_dumper)
else:
	options = dict()


from .common import *
from .transfer import transfer, transfer_dump, transfer_load
from .model import model

def transfer_init_and_dump(transfer_functions):
	transfer_functions[channel_dict['H-Ion']] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch1.dat'))
	transfer_dump(transfer_functions[channel_dict['H-Ion']], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch1.obj'))
	transfer_functions[channel_dict['He-Ion']] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch2.dat'))
	transfer_dump(transfer_functions[channel_dict['He-Ion']], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch2.obj'))
	transfer_functions[channel_dict['Ly-A']] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch3.dat'))
	transfer_dump(transfer_functions[channel_dict['Ly-A']], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch3.obj'))
	transfer_functions[channel_dict['Heat']] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch4.dat'))
	transfer_dump(transfer_functions[channel_dict['Heat']], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch4.obj'))
	transfer_functions[channel_dict['LowE']] = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch5.dat'))
	transfer_dump(transfer_functions[channel_dict['LowE']], os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch5.obj'))

	return transfer_functions

def transfer_load_from_dump(transfer_functions):
	transfer_functions[channel_dict['H-Ion']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch1.obj') )
	transfer_functions[channel_dict['He-Ion']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch2.obj') )
	transfer_functions[channel_dict['Ly-A']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch3.obj') )
	transfer_functions[channel_dict['Heat']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch4.obj') )
	transfer_functions[channel_dict['LowE']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch5.obj') )

	return transfer_functions

#################################

transfer_functions = np.empty(shape=5, dtype=transfer)

transfer_is_initialized = True
for i in xrange(1,6):
	transfer_is_initialized = transfer_is_initialized and os.path.isfile(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch%i.obj' % i))

if not transfer_is_initialized:
	print_info('The transfer seem not to be initialized. This will be done now. this may take a few seconds.')
	transfer_init_and_dump(transfer_functions)
	print_info('The transfer functions are now initialized and loaded.\n')
else:
	#print_info('The transfer functions are already initialized and loaded.\n')
	transfer_load_from_dump(transfer_functions)
del transfer_is_initialized

logEnergies = transfer_functions[0].log10E[:]
redshift = transfer_functions[0].z_deposited[:]
