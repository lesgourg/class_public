# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:25:10 2016

@author: patrick
"""
import sys
import os
import DarkAges
from DarkAges import *

current_dir = os.path.dirname(sys.argv[0])

transfer_functions1 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch1.obj'),'rb') )
transfer_functions2 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch2.obj'),'rb') )
transfer_functions3 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch3.obj'),'rb') )
transfer_functions4 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch4.obj'),'rb') )
transfer_functions5 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch5.obj'),'rb') )

logEnergies = transfer_functions1.log10E[:]

exponent = float(sys.argv[1])
if exponent < 0 or exponent > 3:
	print('BREAK: The mass should be in the range [5 GeV, 5000 GeV]. Hence the exponent needs to be in the range [0, 3].')
	raise SystemExit
sampling_mass = 5*(10**exponent) 
#if len(sys.argv) < 3:
#	annihilation = 1.683e-5 / sampling_mass
#else:
#	annihilation = float(sys.argv[2])*1/sampling_mass

sampling_spectrum_el = np.zeros_like(logEnergies, dtype=np.float64)
sampling_spectrum_ph = np.zeros_like(logEnergies, dtype=np.float64)
sampling_spectrum_oth = np.zeros_like(logEnergies, dtype=np.float64)

spec_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_interp.obj'),'rb')
#spec_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_interp_noEW.obj'),'rb')
#spec_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_interp_EW.obj'),'rb')

## NOTE: 'interp_array' has 3 dimensions (primaries, final states, log10E)
interp_array = dill.load( spec_interp_file )

sampling_spectrum_el = np.zeros_like(logEnergies)
sampling_spectrum_ph = np.zeros_like(logEnergies)
sampling_spectrum_oth = np.zeros_like(logEnergies)

for idx_E in xrange(interp_array.shape[2]):
	temp_el = np.zeros(interp_array.shape[0], dtype=np.float64)
	temp_ph = np.zeros(interp_array.shape[0], dtype=np.float64)
	temp_oth = np.zeros(interp_array.shape[0], dtype=np.float64)
	for idx_primary in xrange(interp_array.shape[0]):
		temp_el[idx_primary] = max(0,interp_array[idx_primary,0,idx_E](sampling_mass))
		temp_ph[idx_primary] = max(0,interp_array[idx_primary,1,idx_E](sampling_mass))
		temp_oth[idx_primary] = max(0,interp_array[idx_primary,2,idx_E](sampling_mass))

	sampling_spectrum_el[idx_E] = sum( temp_el )
	sampling_spectrum_ph[idx_E] = sum( temp_ph )	
	sampling_spectrum_oth[idx_E] = sum( temp_oth )


## Calculate the f-functions for each channel
mix_model = model(sampling_spectrum_el, sampling_spectrum_ph, sampling_spectrum_oth, 1e9*sampling_mass)
# !!Note that the order of the channels differ!!
z, f1 = mix_model.calc_f(transfer_functions1) # Channel: H-Ion
f2 = mix_model.calc_f(transfer_functions2)[-1] # Channel: He-Ion
f3 = mix_model.calc_f(transfer_functions3)[-1] # Channel: Ly-A
f4 = mix_model.calc_f(transfer_functions4)[-1] # Channel: Heating
f5 = mix_model.calc_f(transfer_functions5)[-1] # Channel: LowE

first = 2
last = len(z) - 7
min_z = 0.
max_z = 1e4 

print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\n'%( (last-first) + 2))
#print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\t%.8g\n'%( (last-first) + 2, annihilation))
print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(min_z,f4[first],f3[first],f1[first],f2[first],f5[first]))
for idx in range(first,last):
	print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(z[idx],f4[idx],f3[idx],f1[idx],f2[idx],f5[idx]))
print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(max_z,f4[last-1],f3[last-1],f1[last-1],f2[last-1],f5[last-1]))
#print('\n\nmass: %.2g GeV\n'%sampling_mass)
