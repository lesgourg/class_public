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
from DarkAges import redshift

if len(sys.argv) <3:
	print_error("There are too few arguments passed. I expected at least 2")
sampling_mass = float(sys.argv[1])
if sampling_mass < 50 or sampling_mass > 100:
	print_error("The mass-parameter sholud be in the range [50 GeV, 100 GeV]")
mixing = float(sys.argv[2])
if mixing < 0 or mixing > 1:
	print_error("The mixing-parameter sholud be in the range [0,1]")

model_dir = os.path.split(os.path.realpath(__file__))[0]
model_name =  model_dir.split('/')[-1]

with open(os.path.join(model_dir, '{}.obj'.format(model_name)),'rb') as dump_file:
	dump_dict = dill.load(dump_file)
	interpolated_f = dump_dict.get('f-function')

f_functions = np.zeros(shape=(5,len(redshift)))
for idx_chan, z in itertools.product(*(channel_dict.values(),redshift)):
	idx_z = get_index(redshift, z)
	temp1 = interpolated_f[0,idx_chan,idx_z](sampling_mass) # Primaries are muons
	temp2 = interpolated_f[1,idx_chan,idx_z](sampling_mass) # Primaries are bottom
	f_functions[idx_chan, idx_z] = mixing * temp1 + (1-mixing)*temp2

finalize(redshift, 
         f_functions[channel_dict['Heat']], 
         f_functions[channel_dict['Ly-A']], 
         f_functions[channel_dict['H-Ion']],
         f_functions[channel_dict['He-Ion']],
         f_functions[channel_dict['LowE']])
