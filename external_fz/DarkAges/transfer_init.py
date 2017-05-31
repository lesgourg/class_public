import os
from common import *
from transfer import transfer

file_ch1 = open(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch1.obj'),'wb')
transfer_functions = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch1.dat'))
dill.dump(transfer_functions, file_ch1)
file_ch1.close()
del transfer_functions

file_ch2 = open(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch2.obj'),'wb')
transfer_functions = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch2.dat'))
dill.dump(transfer_functions, file_ch2)
file_ch2.close()
del transfer_functions

file_ch3 = open(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch3.obj'),'wb')
transfer_functions = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch3.dat'))
dill.dump(transfer_functions, file_ch3)
file_ch3.close()
del transfer_functions

file_ch4 = open(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch4.obj'),'wb')
transfer_functions = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch4.dat'))
dill.dump(transfer_functions, file_ch4)
file_ch4.close()
del transfer_functions

file_ch5 = open(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/transfer_Ch5.obj'),'wb')
transfer_functions = transfer(os.path.join(os.environ['DARKAGES_BASE'],'transfer_functions/original/Transfer_Ch5.dat'))
dill.dump(transfer_functions, file_ch5)
file_ch5.close()
del transfer_functions

