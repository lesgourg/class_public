# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:25:10 2016

@author: patrick
"""
import sys
import DarkAges_Class
from DarkAges_Class import *

if (len(sys.argv) != 3):
	print "Method %s takes exactly two arguments. Please provide a file containing the spectrum and a file to save the f-functions and call the method in the following scheme '%s Input.dat Output.dat'"% (sys.argv[0],sys.argv[0])
	raise SystemExit

infile = sys.argv[1]
outfile = sys.argv[2]

print "Running with input: %s"%sys.argv[1]

transfer_functions1 = dill.load( open('./transfer_Ch1.obj') )
transfer_functions2 = dill.load( open('./transfer_Ch2.obj') )
transfer_functions3 = dill.load( open('./transfer_Ch3.obj') )
transfer_functions4 = dill.load( open('./transfer_Ch4.obj') )
transfer_functions5 = dill.load( open('./transfer_Ch5.obj') )

# Muon-spectrum (M_DM = 100 GeV)
spec_data = np.genfromtxt(infile, unpack=True, usecols=(0,1,2,3,4), skip_header=1)

read_log10energy = 9*np.ones_like(spec_data[1,:])+spec_data[1,:]
read_el = spec_data[2,:]
read_ph = spec_data[3,:]
read_oth = spec_data[4,:]

m = spec_data[0,0]*1e9
spectrum_el = log_fit(read_log10energy, 1e-9*read_el, transfer_functions1.log10E)
spectrum_ph = log_fit(read_log10energy, 1e-9*read_ph, transfer_functions1.log10E)
spectrum_oth = log_fit(read_log10energy, 1e-9*read_oth, transfer_functions1.log10E)

## Rescale the interpolated spectra
rescaling_e = trapz( logConversion(transfer_functions1.log10E)*spectrum_el, logConversion(transfer_functions1.log10E) ) /  trapz( logConversion(read_log10energy)*(read_el)*1e-9, logConversion(read_log10energy) )
rescaling_p = trapz( logConversion(transfer_functions1.log10E)*spectrum_ph, logConversion(transfer_functions1.log10E) ) /  trapz( logConversion(read_log10energy)*(read_ph)*1e-9, logConversion(read_log10energy) )
rescaling_o = trapz( logConversion(transfer_functions1.log10E)*spectrum_oth, logConversion(transfer_functions1.log10E) ) /  trapz( logConversion(read_log10energy)*(read_oth)*1e-9, logConversion(read_log10energy) )

spectrum_el *= 1/rescaling_e
spectrum_ph *= 1/rescaling_p
spectrum_oth *= 1/rescaling_o
##

model = model(spectrum_el,spectrum_ph,spectrum_oth,m)
z, f1 = model.calc_f(transfer_functions1) # Channel: H-Ion
f2 = model.calc_f(transfer_functions2)[-1] # Channel: He-Ion
f3 = model.calc_f(transfer_functions3)[-1] # Channel: Ly-A
f4 = model.calc_f(transfer_functions4)[-1] # Channel: Heating
f5 = model.calc_f(transfer_functions5)[-1] # Channel: LowE

last = len(z) - 6

file_out = open(outfile, 'w')
file_out.write('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n%i\n'%(last + 2))
file_out.write('\n%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(0.,f4[0],f3[0],f1[0],f2[0],f5[0]))
for idx in range(last):
	file_out.write('\n%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(z[idx],f4[idx],f3[idx],f1[idx],f2[idx],f5[idx]))
file_out.write('\n%.4g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g'%(10000.,f4[last-1],f3[last-1],f1[last-1],f2[last-1],f5[last-1]))
file_out.close()
#print 'Saved f(z)-curves under "%s"'%outfile

