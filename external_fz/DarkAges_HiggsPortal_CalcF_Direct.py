# -*- coding: utf-8 -*-

## Order of the primaries:
# 0:'top', 1:'higgs', 2:'zboson', 3:'wboson', 4:'bottom', 5:'tau', 6:'charm', 7:'gluon', 8:'gamma'

import sys
import os
import DarkAges
from DarkAges import *

current_dir = os.path.dirname(sys.argv[0])

#mass = float(sys.argv[1])
#if mass < 5e1 or mass > 5e3:
#	print "BREAK: The mass should be in the range [5 GeV, 5000 GeV]."
#	raise SystemExit

log10_of_mass_div5 = float(sys.argv[1])
if log10_of_mass_div5 < 0 or log10_of_mass_div5 > 3:
	print "BREAK: The mass should be in the range [5 GeV, 5000 GeV]. Therefore you have to choose the input parameter (log10 of m / 5) to be in the range [0,3]."
	raise SystemExit
mass = 5 * 10**log10_of_mass_div5

mass = np.asarray(mass)

transfer_functions = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch1.obj'),'rb') )
redshift = transfer_functions.z_deposited[:]
del transfer_functions

fgrid_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_fgrid_interp.obj'),'rb')
frac_interp_file = open(os.path.join(current_dir, 'spectra/HiggsPortal_fractions_interp.obj'),'rb')

## NOTE: 'fgrid_interp_array' has 3 dimensions (primaries, deposition channel, redshift)
## NOTE: 'frac_interp_array' has 1 dimension (primaries)
fgrid_interp_array = dill.load( fgrid_interp_file )
frac_interp_array = dill.load( frac_interp_file )

fgrid_interp_file.close()
frac_interp_file.close()

assert frac_interp_array.shape[0] == fgrid_interp_array.shape[0]

f_function = np.zeros(shape=(fgrid_interp_array.shape[1], fgrid_interp_array.shape[2]), dtype=np.float64)
for idx_Z in xrange(fgrid_interp_array.shape[2]):
	for idx_ch in xrange(fgrid_interp_array.shape[1]):
		temp_f_at_z = np.zeros(shape=(fgrid_interp_array.shape[0]), dtype=np.float64)
		fractions = np.zeros(shape=(fgrid_interp_array.shape[0]), dtype=np.float64)
		for idx_primary in xrange(fgrid_interp_array.shape[0]):
			temp_f_at_z[idx_primary] = fgrid_interp_array[idx_primary,idx_ch,idx_Z](mass)
			fractions[idx_primary] = frac_interp_array[idx_primary](mass)
		frac_norm = sum(fractions)
		#print frac_norm 
		if frac_norm > 0:
			fractions *= 1/frac_norm 
		f_function[idx_ch,idx_Z] = np.inner( temp_f_at_z, fractions )

f0 = f_function[0,:]
f1 = f_function[1,:]
f2 = f_function[2,:]
f3 = f_function[3,:]
f4 = f_function[4,:]
f5 = f_function[5,:]

first = 2
last = len(redshift) - 7
min_z = 0.
max_z = 1e4 

#print('#z_dep\tf_tot\tf_1t\tf_2\tf_3\tf_4\tf_5\n\n%i\n'%( (last-first) + 2))
#print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(min_z,f0[first],f1[first],f2[first],f3[first],f4[first],f5[first]))
#for idx in range(first,last):
#	print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(redshift[idx],f0[idx],f1[idx],f2[idx],f3[idx],f4[idx],f5[idx]))
#print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(max_z,f0[last-1],f1[last-1],f2[last-1],f3[last-1],f4[last-1],f5[last-1]))

print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\n'%( (last-first) + 2))
#print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\t%.8g\n'%( (last-first) + 2, annihilation))
print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(min_z,f4[first],f3[first],f1[first],f2[first],f5[first]))
for idx in range(first,last):
	print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(redshift[idx],f4[idx],f3[idx],f1[idx],f2[idx],f5[idx]))
print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(max_z,f4[last-1],f3[last-1],f1[last-1],f2[last-1],f5[last-1]))


