import os
from common import *
from transfer import transfer, transfer_load

transfer_functions = np.empty(shape=5, dtype=transfer)

transfer_functions[channel_dict['H-Ion']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch1.obj') )	
transfer_functions[channel_dict['He-Ion']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch2.obj') )	
transfer_functions[channel_dict['Ly-A']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch3.obj') )	
transfer_functions[channel_dict['Heat']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch4.obj') )	
transfer_functions[channel_dict['LowE']] = transfer_load( os.path.join(os.environ['DARKAGES_BASE'], 'transfer_functions/transfer_Ch5.obj') )	

logEnergies = transfer_functions[0].log10E[:]
redshift = transfer_functions[0].z_deposited[:]
