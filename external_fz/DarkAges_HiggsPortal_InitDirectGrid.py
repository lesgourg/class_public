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
	mass_mask = ( abs(spec_data[0,:]/float(mass) - 1) < 1e-5 ) 

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

	## METHOD 1: Interpolate the dN/dE spectrum
	spectrum_el = log_fit(read_log10energy, 1e-9*read_el/rescaling, sampling_log10E)
	spectrum_ph = log_fit(read_log10energy, 1e-9*read_ph/rescaling, sampling_log10E)
	spectrum_oth = log_fit(read_log10energy, 1e-9*read_oth/rescaling, sampling_log10E)

	### METHOD 2: Interpolate the E.dN/dE spectrum
	#spectrum_el = log_fit(read_log10energy, 1e-9*read_el*logConversion(read_log10energy)/rescaling, sampling_log10E)
	#spectrum_ph = log_fit(read_log10energy, 1e-9*read_ph*logConversion(read_log10energy)/rescaling, sampling_log10E)
	#spectrum_oth = log_fit(read_log10energy, 1e-9*read_oth*logConversion(read_log10energy)/rescaling, sampling_log10E)

	#spectrum_el *= 1/logConversion(sampling_log10E)
	#spectrum_ph *= 1/logConversion(sampling_log10E)
	#spectrum_oth *= 1/logConversion(sampling_log10E)

	## Rescale the interpolated spectra
	if non_zero_spec:
		rescaling_e = trapz( logConversion(sampling_log10E)*spectrum_el, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_el/rescaling)*1e-9, logConversion(read_log10energy) )
		rescaling_p = trapz( logConversion(sampling_log10E)*spectrum_ph, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_ph/rescaling)*1e-9, logConversion(read_log10energy) )
		rescaling_o = trapz( logConversion(sampling_log10E)*spectrum_oth, logConversion(sampling_log10E) ) /  trapz( logConversion(read_log10energy)*(read_oth/rescaling)*1e-9, logConversion(read_log10energy) )

		spectrum_el *= 1/rescaling_e
		spectrum_ph *= 1/rescaling_p
		spectrum_oth *= 1/rescaling_o

	for idx in xrange(len(spectrum_el)):
		spectrum_el[idx] = max(0, spectrum_el[idx])
		spectrum_ph[idx] = max(0, spectrum_ph[idx])
		spectrum_oth[idx] = max(0, spectrum_oth[idx])

	return np.array([spectrum_el, spectrum_ph, spectrum_oth])

def getF(spec_el, spec_ph, spec_oth, trans_instance, mass):
	tmp_model = model(spec_el, spec_ph, spec_oth, 1e9*mass)
	f_function = tmp_model.calc_f(trans_instance)[-1]
	
	first = 2
	last = len(f_function) - 7

	ret_f = np.zeros_like(f_function)
	
	ret_f[:first] = f_function[first]
	ret_f[first:last] = f_function[first:last]
	ret_f[last:] = f_function[last]

	return ret_f

#######################
#######################
#######################

current_dir = os.path.dirname(sys.argv[0])
#print current_dir

transfer_functions1 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch1.obj'),'rb') )
transfer_functions2 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch2.obj'),'rb') )
transfer_functions3 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch3.obj'),'rb') )
transfer_functions4 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch4.obj'),'rb') )
transfer_functions5 = dill.load( open(os.path.join(current_dir, 'transfer_functions/transfer_Ch5.obj'),'rb') )

logEnergies = transfer_functions1.log10E[:]
redshift = transfer_functions1.z_deposited[:]

masses = np.genfromtxt(os.path.join(current_dir, 'spectra/Cirelli_masses.dat'),unpack=True, usecols=(0), dtype=np.float64)
masses = masses[masses <= 5000]

#primaries = np.array(['top', 'higgs', 'zzstar_merged', 'wwstar_merged', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a16'))
primaries = np.array(['top', 'higgs', 'zzstar_merged2', 'wwstar_merged2', 'bottom', 'tau', 'charm', 'gluon', 'gamma'], dtype=np.dtype('a16'))


specfiles = np.empty(shape=(len(primaries)), dtype=np.dtype('a128'))
for idx, primary in enumerate(primaries):
	specfiles[idx] = os.path.join(current_dir,'spectra/%s_EW.dat' % primary) 

frac_read = np.genfromtxt(os.path.join(current_dir,'spectra/HiggsPortal_fractions.dat'), unpack=True, usecols=(0,4,5,6,7,8,9,10,11,12), dtype=np.float64)
fractions = np.empty(shape=(len(primaries),len(masses)), dtype=np.float64)
fractions_interp_obj = np.empty(shape=(len(primaries)), dtype=logLinearInterpolator)
#fractions_interp_obj = np.empty(shape=(len(primaries)), dtype=logInterpolator)
#fractions_interp_obj = np.empty(shape=(len(primaries)), dtype=interp1d)
for idx in xrange(len(fractions[:,0])):
	temp_interp = logLinearInterpolator(frac_read[0,:], frac_read[idx+1,:], 0)
	#temp_interp = logInterpolator(frac_read[0,:], frac_read[idx+1,:], 0)
	#temp_interp = interp1d(frac_read[0,:], frac_read[idx+1,:], kind='linear')
	fractions[idx,:] = temp_interp(masses)
	fractions_interp_obj[idx] = temp_interp
	#print fractions[idx,:]
del temp_interp

frac_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_fractions_interp.obj'),'wb')
dill.dump(fractions_interp_obj, frac_interp_file)
frac_interp_file.close()

# Reweighting the fractions
for idx in xrange(fractions.shape[-1]):
	reweight = sum(fractions[:,idx])
	fractions[:,idx] *= 1/reweight

print "Interpolating the fractions. DONE"

progress = 0.
num_of_steps = len(masses)

f_grid = np.zeros( shape=(len(masses),len(primaries),len(redshift),6), dtype=np.float64 )
for idx_mass, mass in enumerate(masses):
	progress += 100./(num_of_steps)
	for idx_primary, fname in enumerate(specfiles):
		tmp_spec =  np.zeros( shape=(3,len(logEnergies)), dtype=np.float64 )
		tmp_spec[:,:] = getSpec(fname, mass, logEnergies)
		f_grid[idx_mass,idx_primary,:,1] = getF(tmp_spec[0,:], tmp_spec[1,:], tmp_spec[2,:], transfer_functions1, mass)
		f_grid[idx_mass,idx_primary,:,2] = getF(tmp_spec[0,:], tmp_spec[1,:], tmp_spec[2,:], transfer_functions2, mass)
		f_grid[idx_mass,idx_primary,:,3] = getF(tmp_spec[0,:], tmp_spec[1,:], tmp_spec[2,:], transfer_functions3, mass)
		f_grid[idx_mass,idx_primary,:,4] = getF(tmp_spec[0,:], tmp_spec[1,:], tmp_spec[2,:], transfer_functions4, mass)
		f_grid[idx_mass,idx_primary,:,5] = getF(tmp_spec[0,:], tmp_spec[1,:], tmp_spec[2,:], transfer_functions5, mass)
		f_grid[idx_mass,idx_primary,:,0] = np.sum(f_grid[idx_mass,idx_primary,:,1:], axis=1)
	print '\r(Progess: %.2f %%)' % (progress)

print "Done: Preparing the f-function grid"


#####
# Interpolate the f-function grid.

#interp_array = np.zeros( shape=(len(primaries),6,len(redshift)) ,dtype=logLinearInterpolator )
#for idx_Z in xrange(len(redshift)):
#	for idx_ch in xrange(6):
#		for idx_primary in xrange(len(primaries)):
#			#FIXME
#			interp_array[idx_primary,idx_ch,idx_Z] = logLinearInterpolator(masses, f_grid[:,idx_primary,idx_Z,idx_ch], 0)
#			#FIXME

interp_array = np.zeros( shape=(len(primaries),6,len(redshift)) ,dtype=logInterpolator )
for idx_Z in xrange(len(redshift)):
	for idx_ch in xrange(6):
		for idx_primary in xrange(len(primaries)):
			#FIXME
			interp_array[idx_primary,idx_ch,idx_Z] = logInterpolator(masses, f_grid[:,idx_primary,idx_Z,idx_ch], 0)
			#FIXME

#Dill the data away:
fgrid_interp_file =  open(os.path.join(current_dir, 'spectra/HiggsPortal_fgrid_interp.obj'),'wb')
dill.dump(interp_array, fgrid_interp_file)
fgrid_interp_file.close()

del interp_array
print "Interpolation ( Sliced 1D ) done!"


