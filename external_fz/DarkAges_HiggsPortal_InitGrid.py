# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:25:10 2016

@author: patrick
"""
import sys
import os
import DarkAges
from DarkAges import *

def getSpec(filename, mass, sampling_log10E):
	spec_data = np.genfromtxt(filename, unpack=True, usecols=(0,5,6,7,8), skip_header=1, dtype=np.float64)
	mass_mask = ( abs(spec_data[0,:]/float(mass) - 1) < 1e-3 ) 
	if len(spec_data[0,mass_mask]) == 0:
		print "something went wrong with the mass_mask"

	read_log10energy = 9*np.ones_like(spec_data[1,mass_mask])+spec_data[1,mass_mask]
	read_el = spec_data[2,mass_mask]
	read_ph = spec_data[3,mass_mask]
	read_oth = spec_data[4,mass_mask]

	m = mass*1e9
	non_zero_spec = True
	rescaling = trapz( logConversion(read_log10energy)*(read_el+read_ph+read_oth)*1e-9, logConversion(read_log10energy) ) / (2*m)
	if rescaling < 1e-6:
		rescaling = 1.
		non_zero_spec = False
	#rescaling = 1

	spectrum_el = log_fit(read_log10energy, 1e-9*read_el/rescaling, sampling_log10E)
	spectrum_ph = log_fit(read_log10energy, 1e-9*read_ph/rescaling, sampling_log10E)
	spectrum_oth = log_fit(read_log10energy, 1e-9*read_oth/rescaling, sampling_log10E)

	# In the unphysical region E > m the spectrum needs to vanish
	unphysical_region_mask =  (logConversion(sampling_log10E) > m)
	spectrum_el[unphysical_region_mask] = 0.
	spectrum_ph[unphysical_region_mask] = 0.
	spectrum_oth[unphysical_region_mask] = 0.

	for idx in xrange(len(spectrum_el)):
		spectrum_el[idx] = max(0., spectrum_el[idx])
		spectrum_ph[idx] = max(0., spectrum_ph[idx])
		spectrum_oth[idx] = max(0., spectrum_oth[idx])	

	## Rescale the interpolated spectra
	if non_zero_spec:
		rescaling_e = trapz( logConversion(sampling_log10E)*spectrum_el, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_el/rescaling)*1e-9, logConversion(read_log10energy) )
		rescaling_p = trapz( logConversion(sampling_log10E)*spectrum_ph, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_ph/rescaling)*1e-9, logConversion(read_log10energy) )
		rescaling_o = trapz( logConversion(sampling_log10E)*spectrum_oth, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_oth/rescaling)*1e-9, logConversion(read_log10energy) )

		spectrum_el *= 1/rescaling_e
		spectrum_ph *= 1/rescaling_p
		spectrum_oth *= 1/rescaling_o

	return np.array([spectrum_el, spectrum_ph, spectrum_oth])

#######################
#######################
#######################

current_dir = os.path.dirname(sys.argv[0])

transfer_functions = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch1.obj'),'rb') )
logEnergies = transfer_functions.log10E[:]
del transfer_functions

masses = np.genfromtxt(os.path.join(current_dir, 'spectra/Cirelli_masses.dat'),unpack=True, usecols=(0), dtype=np.float64)
masses = masses[masses <= 5000]

#primaries = np.array(['top', 'higgs', 'zzstar_1', 'wwstar_1', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a16'))
#primaries = np.array(['top', 'higgs', 'zzstar', 'wwstar', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a16'))
#primaries = np.array(['top', 'higgs', 'zzstar_merged', 'wwstar_merged', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a16'))
primaries = np.array(['top', 'higgs', 'zzstar_merged2', 'wwstar_merged2', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a16'))

specfiles = np.empty(len(primaries), dtype=np.dtype('a128'))
for idx, primary in enumerate(primaries):
	specfiles[idx] = os.path.join(current_dir, 'spectra/%s_noEW.dat' % primary) 

frac_read = np.genfromtxt(os.path.join(current_dir, 'spectra/HiggsPortal_fractions.dat'), unpack=True, usecols=(0,4,5,6,7,8,9,10,11,12), dtype=np.float64)
fractions = np.empty(shape=(len(primaries),len(masses)), dtype=np.float64)
for idx in xrange(len(fractions[:,0])):
	temp_interp = interp1d(frac_read[0,:], frac_read[idx+1,:], kind='linear')
	fractions[idx,:] = temp_interp(masses)
del temp_interp

# Reweighting the fractions
#for idx in xrange(fractions.shape[-1]):
#	reweight = sum(fractions[:,idx])
#	fractions[:,idx] *= 1/reweight

print "Done: Interpolating the fractions."

spectra_grid = np.zeros( shape=(len(masses),len(primaries),3,len(logEnergies)), dtype=np.float64 )
for idx_mass, mass in enumerate(masses):
	for idx_primary, fname in enumerate(specfiles):
		spectra_grid[idx_mass,idx_primary,:,:] = getSpec(fname, mass, logEnergies) 

print "Done: Preparing the Spectra"
#####
# Interpolate the noEW spectra.
interp_array = np.zeros( shape=(len(primaries),3,len(logEnergies)) ,dtype=interp1d)
for idx_E in xrange(len(logEnergies)):
	for idx_primary in xrange(len(primaries)):
		interp_array[idx_primary,0,idx_E] = interp1d(masses, fractions[idx_primary,:]*spectra_grid[:,idx_primary,0,idx_E], kind='quadratic')
		interp_array[idx_primary,1,idx_E] = interp1d(masses, fractions[idx_primary,:]*spectra_grid[:,idx_primary,1,idx_E], kind='quadratic')
		interp_array[idx_primary,2,idx_E] = interp1d(masses, fractions[idx_primary,:]*spectra_grid[:,idx_primary,2,idx_E], kind='quadratic')
#Dill the data away:
spec_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_interp.obj'),'wb')
dill.dump(interp_array, spec_interp_file)
spec_interp_file.close()
del interp_array
print "Interpolation done."

