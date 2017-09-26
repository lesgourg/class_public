u"""
.. module:: DarkAges
   :synopsis: Module for the calculation of f(z) given a certain spectrum of particles
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

When loaded by :code:`import DarkAges` or by :code:`from DarkAges import ...`
"""

import numpy as np
np.seterr(all='ignore')
#Silence numpy when it comes to overflows etc.
import sys
import os

class DarkAgesError(Exception):
	def __init__(self, message):
		Exception.__init__(self)
		self.message = '\n\n!!! ERROR !!!\n\n --> {} \n'.format(message)
		#print('\n\n ERROR:\n !!! {} !!!\n'.format(message))

	def __str__(self):
		return self.message

os.environ['DARKAGES_BASE'] = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.insert(0, os.environ['DARKAGES_BASE'])

DarkOptions = dict()
DarkOptions.update({'background_H0':67.27, 'background_Omega_m': 0.3156, 'background_Omega_r' :8e-5})
if 'DARKAGES_TOPLEVEL_PID' in os.environ:
	if os.getppid() == int(os.environ['DARKAGES_TOPLEVEL_PID']):
		import yaml
		top_level_random = int(os.environ['DARKAGES_TOPLEVEL_RANDN'])
		with open(os.path.join( os.environ['DARKAGES_BASE'], 'pid_{:d}_{:d}.yaml'.format(os.getppid(), top_level_random)), 'rb') as options_dumper:
			loaded_options = yaml.load(options_dumper)
			DarkOptions.update( loaded_options )
			del loaded_options
			del top_level_random

from .common import print_info, channel_dict
from .transfer import transfer
#from .model import model

def _transfer_init_and_dump(transfer_functions):
	from .transfer import transfer_dump
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

def _transfer_load_from_dump(transfer_functions):
	from .transfer import transfer_load
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
	_transfer_init_and_dump(transfer_functions)
	print_info('The transfer functions are now initialized and loaded.\n')
else:
	#print_info('The transfer functions are already initialized and loaded.\n')
	_transfer_load_from_dump(transfer_functions)
# def set_energy_range(Emin,Emax,nbins):
#     Emin_table = Emin
#     Emax_table = Emax
#     nbins_table = nbins
# Emin_table = 0.1
# Emax_table = 12.7675
Emin_table = None
Emax_table = None
nbins_table = 100
# LogEnergies = np.logspace(Emin_table,Emax_table,nbins_table)
if Emin_table == None or Emax_table == None:
    logEnergies = transfer_functions[0].log10E[:]
else:
    logEnergies = np.linspace(Emin_table,Emax_table,nbins_table)
# print logEnergies
redshift = transfer_functions[0].z_deposited[:]

# Clean from functions and variables which should not be available when
# "DarkAges" is imported.
del transfer_is_initialized
del i
del print_info
