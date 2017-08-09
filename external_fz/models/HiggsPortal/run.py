'''
This model has two input parameters:
    1) mass (between 50 and 100 GeV)
    2) the relative contribution of muons as primaries	 
'''

import numpy as np
import itertools
import dill
import os 
import sys
if os.environ['DARKAGES_BASE']:
	sys.path.insert(0, os.environ['DARKAGES_BASE'] )
import DarkAges
from DarkAges.common import print_error, finalize, channel_dict, get_index
from DarkAges import redshift, options

if len(sys.argv) <2:
	print_error("There are too few arguments passed. I expected at least 1")
sampling_mass = float(sys.argv[1])
if sampling_mass < 5 or sampling_mass > 5e3:
	print_error("The mass-parameter sholud be in the range [5 GeV, 5 TeV]")

model_dir = os.path.split(os.path.realpath(__file__))[0]
model_name =  model_dir.split('/')[-1]

with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'rb') as dump_file:
	dump_dict = dill.load(dump_file)
	interpolated_f = dump_dict.get('f-functions')
	interpolated_frac = dump_dict.get('fractions')

f_functions = np.zeros(shape=(5,len(redshift)))
for idx_chan, z in itertools.product(*(channel_dict.values(),redshift)):
	idx_z = get_index(redshift, z)
	temp_f = np.zeros(shape=(interpolated_f.shape[0]), dtype=np.float64)
	temp_frac = np.zeros(shape=(interpolated_f.shape[0]), dtype=np.float64)
	for idx_prim in xrange(interpolated_f.shape[0]):
		temp_f[idx_prim] = interpolated_f[idx_prim,idx_chan,idx_z](sampling_mass)
		temp_frac[idx_prim] = interpolated_frac[idx_prim](sampling_mass)
	f_functions[idx_chan, idx_z] = np.tensordot(temp_f, temp_frac, axes=(0,0))

finalize(redshift, 
         f_functions[channel_dict['Heat']], 
         f_functions[channel_dict['Ly-A']], 
         f_functions[channel_dict['H-Ion']],
         f_functions[channel_dict['He-Ion']],
         f_functions[channel_dict['LowE']])
