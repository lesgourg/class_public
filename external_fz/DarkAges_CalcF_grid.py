# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:25:10 2016

@author: patrick
"""
import sys
import os
import DarkAges
from DarkAges import *

def getSpec(filename, sampling_log10E):
	spec_data = np.genfromtxt(filename, unpack=True, usecols=(0,1,2,3,4), skip_header=1)

	read_log10energy = 9*np.ones_like(spec_data[1,:])+spec_data[1,:]
	read_el = spec_data[2,:]
	read_ph = spec_data[3,:]
	read_oth = spec_data[4,:]

	m = spec_data[0,0]*1e9
	rescaling = trapz( logConversion(read_log10energy)*(read_el+read_ph+read_oth)*1e-9, logConversion(read_log10energy) ) / (2*m)
	#rescaling = 1

	spectrum_el = log_fit(read_log10energy, 1e-9*read_el/rescaling, sampling_log10E)
	spectrum_ph = log_fit(read_log10energy, 1e-9*read_ph/rescaling, sampling_log10E)
	spectrum_oth = log_fit(read_log10energy, 1e-9*read_oth/rescaling, sampling_log10E)

	## Rescale the interpolated spectra
	rescaling_e = trapz( logConversion(sampling_log10E)*spectrum_el, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_el)*1e-9, logConversion(read_log10energy) )
	rescaling_p = trapz( logConversion(sampling_log10E)*spectrum_ph, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_ph)*1e-9, logConversion(read_log10energy) )
	rescaling_o = trapz( logConversion(sampling_log10E)*spectrum_oth, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_oth)*1e-9, logConversion(read_log10energy) )

	spectrum_el *= 1/rescaling_e
	spectrum_ph *= 1/rescaling_p
	spectrum_oth *= 1/rescaling_o

	return np.array([spectrum_el, spectrum_ph, spectrum_oth])

#######################
#######################
#######################

current_dir = os.path.dirname(sys.argv[0])
sampling_mass = float(sys.argv[1])
mixing = float(sys.argv[2])
if mixing < 0 or mixing > 1:
	print "BREAK: The mixing-parameter sholud be in the range [0,1]"
	raise SystemExit
#outfile = sys.argv[3]

trans_file1 =  open(os.path.join(current_dir, 'transfer_functions/transfer_Ch1.obj'))
transfer_functions1 = dill.load( trans_file1 )
trans_file2 =  open(os.path.join(current_dir, 'transfer_functions/transfer_Ch2.obj'))
transfer_functions2 = dill.load( trans_file2 )
trans_file3 =  open(os.path.join(current_dir, 'transfer_functions/transfer_Ch3.obj'))
transfer_functions3 = dill.load( trans_file3 )
trans_file4 =  open(os.path.join(current_dir, 'transfer_functions/transfer_Ch4.obj'))
transfer_functions4 = dill.load( trans_file4 )
trans_file5 =  open(os.path.join(current_dir, 'transfer_functions/transfer_Ch5.obj'))
transfer_functions5 = dill.load( trans_file5 )

logEnergies = transfer_functions1.log10E[:]

#masses = np.array([50,55,60,65,70,80,85,90,95,100]) # exclude 75 for testing quality of interpolation
#masses = np.array([50,55,60,65,70,75,80,85,90,95,100]) # full
masses = np.array([50,60,70,80,90,100])
spectra_grid1 = np.zeros( shape=(len(masses),3,len(logEnergies)) )
spectra_grid2 = np.zeros( shape=(len(masses),3,len(logEnergies)) )
for idx, mass in enumerate(masses):
	fname1 = os.path.join(current_dir,"spectra/muon_%i_spec.dat" % mass)
	fname2 = os.path.join(current_dir,"spectra/bottom_%i_spec.dat" % mass)
	spectra_grid1[idx,:,:] = getSpec(fname1, logEnergies) 
	spectra_grid2[idx,:,:] = getSpec(fname2, logEnergies)

sampling_spectrum_el = np.zeros_like(transfer_functions1.log10E)
sampling_spectrum_ph = np.zeros_like(transfer_functions1.log10E)
sampling_spectrum_oth = np.zeros_like(transfer_functions1.log10E)

for idx in xrange(len(logEnergies)):
	interp_el_m= interp1d(masses, spectra_grid1[:,0,idx], kind='quadratic')
	interp_el_b = interp1d(masses, spectra_grid2[:,0,idx], kind='quadratic')
	interp_ph_m = interp1d(masses, spectra_grid1[:,1,idx], kind='quadratic')
	interp_ph_b = interp1d(masses, spectra_grid2[:,1,idx], kind='quadratic')
	interp_oth_m = interp1d(masses, spectra_grid1[:,2,idx], kind='quadratic')
	interp_oth_b = interp1d(masses, spectra_grid2[:,2,idx], kind='quadratic')
	sampling_spectrum_el[idx] =  mixing*interp_el_m(sampling_mass) + (1 -mixing)*interp_el_b(sampling_mass)
	sampling_spectrum_ph[idx] =  mixing*interp_ph_m(sampling_mass) + (1 -mixing)*interp_ph_b(sampling_mass)	
	sampling_spectrum_oth[idx] =  mixing*interp_oth_m(sampling_mass) + (1 -mixing)*interp_oth_b(sampling_mass)
del interp_el_m, interp_el_b, interp_ph_m, interp_ph_b, interp_oth_m, interp_oth_b
del spectra_grid1, spectra_grid2

## Calculate the f-functions for each channel
mix_model = model(sampling_spectrum_el, sampling_spectrum_ph, sampling_spectrum_oth, 1e9*sampling_mass)
# !!Note that the order of the channels differ!!
# Channel: H-Ion
z, f1 = mix_model.calc_f(transfer_functions1) 
# Channel: He-Ion
f2 = mix_model.calc_f(transfer_functions2)[-1] 
# Channel: Ly-A
f3 = mix_model.calc_f(transfer_functions3)[-1] 
# Channel: Heating
f4 = mix_model.calc_f(transfer_functions4)[-1] 
# Channel: LowE
f5 = mix_model.calc_f(transfer_functions5)[-1] 


first = 2
last = len(z) - 7

#file_out = open(outfile, 'w')
#file_out.write('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\n'%( (last-first) + 2))
#file_out.write('\n%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(0.,f4[first],f3[first],f1[first],f2[first],f5[first]))
#for idx in range(last):
#	file_out.write('\n%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(z[idx],f4[idx],f3[idx],f1[idx],f2[idx],f5[idx]))
#file_out.write('\n%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(10000.,f4[last-1],f3[last-1],f1[last-1],f2[last-1],f5[last-1]))
#file_out.close()
#print 'Saved f(z)-curves under "%s"'%outfile

print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\n'%( (last-first) + 2))
print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(0.,f4[first],f3[first],f1[first],f2[first],f5[first]))
for idx in range(first,last):
	print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(z[idx],f4[idx],f3[idx],f1[idx],f2[idx],f5[idx]))
print('%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(10000.,f4[last-1],f3[last-1],f1[last-1],f2[last-1],f5[last-1]))
