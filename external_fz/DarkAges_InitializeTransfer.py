# -*- coding: utf-8 -*-

import DarkAges
import os
import sys
from DarkAges import *

module_dir = os.path.dirname(sys.argv[0])

#file_all = open(os.path.join(module_dir,'transfer_functions/transfer_all.obj'),'wb')
#transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_all.dat'))
#transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_sum.dat'))
#dill.dump(transfer_functions, file_all)
#file_all.close()
#del transfer_functions

file_ch1 = open(os.path.join(module_dir,'transfer_functions/transfer_Ch1.obj'),'wb')
transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_Ch1.dat'))
dill.dump(transfer_functions, file_ch1)
file_ch1.close()
del transfer_functions

file_ch2 = open(os.path.join(module_dir,'transfer_functions/transfer_Ch2.obj'),'wb')
transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_Ch2.dat'))
dill.dump(transfer_functions, file_ch2)
file_ch2.close()
del transfer_functions

file_ch3 = open(os.path.join(module_dir,'transfer_functions/transfer_Ch3.obj'),'wb')
transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_Ch3.dat'))
dill.dump(transfer_functions, file_ch3)
file_ch3.close()
del transfer_functions

file_ch4 = open(os.path.join(module_dir,'transfer_functions/transfer_Ch4.obj'),'wb')
transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_Ch4.dat'))
dill.dump(transfer_functions, file_ch4)
file_ch4.close()
del transfer_functions

file_ch5 = open(os.path.join(module_dir,'transfer_functions/transfer_Ch5.obj'),'wb')
transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_Ch5.dat'))
dill.dump(transfer_functions, file_ch5)
file_ch5.close()
del transfer_functions

#file_all_fine = open(os.path.join(module_dir,'transfer_functions/transfer_all_fine.obj'),'wb')
#transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_all_fine.dat'))
#dill.dump(transfer_functions, file_all_fine)
#file_all_fine.close()
#del transfer_functions

#file_all_corr1 = open(os.path.join(module_dir,'transfer_functions/transfer_all_corr1.obj'),'wb')
#transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_all_corr1.dat'))
#dill.dump(transfer_functions, file_all_corr1)
#file_all_corr1.close()
#del transfer_functions

#file_all_corr2 = open(os.path.join(module_dir,'transfer_functions/transfer_all_corr2.obj'),'wb')
#transfer_functions = transfer(os.path.join(module_dir,'transfer_functions/original/Transfer_all_corr2.dat'))
#dill.dump(transfer_functions, file_all_corr2)
#file_all_corr2.close()
#del transfer_functions

