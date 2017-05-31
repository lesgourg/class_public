import os
os.environ['DARKAGES_BASE'] = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]

import common
from common import print_info 
from model import model
from transfer import transfer
from interpolator import logInterpolator, logLinearInterpolator

transfer_is_initialized = True
for i in xrange(1,6):
	transfer_is_initialized = transfer_is_initialized and os.path.isfile(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch%i.obj' % i)) 

if not transfer_is_initialized:
	print_info('The transfer seem not to be initialized. This will be done now. this may take a few seconds.')
	from transfer_init import *
	from transfer_load import *
	print_info('The transfer functions are now initialized and loaded.\n')
else:
	#print_info('The transfer functions are already initialized and loaded.\n')
	from transfer_load import *

