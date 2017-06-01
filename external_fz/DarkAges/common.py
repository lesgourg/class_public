from scipy.integrate import trapz
import os
import sys
import numpy as np

###############################
### Definition of functions ###
###############################

''' 
NOTE: the order of the channel given by T. Slatyer and parsed to CLASS differ:
# Slatyer order:  H-Ion, He-Ion, Ly-A,  Heat,   LowE   ->  0,1,2,3,4
# CLASS order:    Heat,  Ly-Aa,  H-Ion, He-Ion, lowE   ->  3,2,0,1,4
'''

channel_dict = {
	  'H-Ion': 0,
	 'He-Ion': 1,
	   'Ly-A': 2,
	   'Heat': 3,
	   'LowE': 4
}

def set_cosmology(H0 = 67.27, Om_M = 0.315, Om_R = 8e-5):
	km_per_Mpc = 3.241e-20

	Cosmo_H0 = H0 * km_per_Mpc # in 1/s
	Cosmo_Omega_M = Om_M
	Cosmo_Omega_R = Om_R
	return Cosmo_H0, Cosmo_Omega_M, Cosmo_Omega_R	

if ('DARKAGES_BG_H0' in os.environ) and ('DARKAGES_BG_OM_M' in os.environ) and ('DARKAGES_BG_OM_R' in os.environ):
	Cosmo_H0, Cosmo_Omega_M, Cosmo_Omega_R = set_cosmology(H0 = float(os.environ['DARKAGES_BG_H0']), 
                                                           Om_M = float(os.environ['DARKAGES_BG_OM_M']), 
                                                           Om_R = float(os.environ['DARKAGES_BG_OM_R']))
else:
	Cosmo_H0, Cosmo_Omega_M, Cosmo_Omega_R = set_cosmology()

def print_info(message):
	print('#INFO: {}'.format(message))

def print_warning(message):
	print('#WARNING: {}'.format(message))

def print_error(message):
	print('\n\n!!! EROOR: {} !!!\n\n'.format(message))
	sys.exit()

def logConversion(LogNumber):
    return 10**LogNumber
    
def simple_H(redshift, H0 = Cosmo_H0, Omega_M = Cosmo_Omega_M):
    '''
    Redshift dependence of H (assuming only matter and radiation to be present). 
    '''
    return H(redshift, Omega_R = 0.)

def H(redshift, H0=Cosmo_H0, Omega_M = Cosmo_Omega_M, Omega_R = Cosmo_Omega_R):
    '''
    Redshift dependence of H (assuming only matter and radiation to be present). 
    '''
    #return H0 * np.sqrt( redshift**(3.) * Omega_M + redshift**(4.) * Omega_R + (1-Omega_R-Omega_M) )
    return H0 * np.sqrt( redshift**(3.) * Omega_M + redshift**(4.) * Omega_R ) * 1/np.sqrt(Omega_R + Omega_M)

def time_at_z(redshift, H0=Cosmo_H0, Omega_M = Cosmo_Omega_M, Omega_R = Cosmo_Omega_R):
    '''
    Time at a given redshift. For this we assume
    only matter and radiation to be present.
    (Taken from EnergyAbsorptionCalculator.nb provided
    as a supplement of arXiV:1211.0283)
    '''
    return 2/(3 * Omega_M**2 * redshift * H0) * ( Omega_M * np.sqrt(Omega_R + (Omega_M / redshift)) + 2 * Omega_R**1.5 * redshift - 2 * Omega_R * np.sqrt(Omega_M + redshift*Omega_R) )

def conversion( redshift, alpha=3 ):
    '''
    This is the conversion factor in front of the z 
    and energy integrals which takes the expnasion into account.
    The alpha is an exponent which depend on the particular case you consider. 
    (For DM-annihilation, this is 3. For DM-Decay it is 0)
    '''
    return ( redshift**alpha / H(redshift) )
    
def nan_clean( array ):
    '''    
    Replacing all NAN-entries of a list with zero. 
    (For example if we sample the injected spectrum 
    for an energy for which we do not have any transfer function)
    '''
    ret = np.empty_like( array, dtype=np.float64 )
    for idx, value in enumerate(array):
        if value == value:
            ret[idx] = value
        else:
            ret[idx] = 0
    return ret

def get_index( array, entry ):
	return np.where(array == entry)[0][0].astype(np.int32)

def f_function_annihilation(logE, z_inj, z_dep, mass_dm, transfer_phot, transfer_elec, 
               spec_phot, spec_elec, spec_tot):
    '''
    This is to calculate the f-function given the photon and electron spectrum
    of the particular model, assuming redshift independent DM-annihalation
    i.e. we assume the thermally averaged crossection to be constant. That is
    why we choose alpha=3 as exponent in the conversion factor.   
    
    The inputs are as follows:
    logE:  
        40x1-dim. object containing log10 of the particle kinetic energies 
        (given in units of eV).
    z_inj / z_dep:
        63x1-dim. object containing the log-spaced redshift (z+1) of 
        deposition (or redshift of injection) 
    transfer_phot / transfer_elec:  
        63x40x63-dim. object containing the transfer function for photons 
        (electrons), taking redshift of deposition, Energy at injection and 
        Redshift of injection. 
    spec_phot / spec_elec:  
        Spectrum (dN/dE) of the injected photons as function of Energy 
        (z dependence already factored out).
    spec_tot:
		Total spectrum (dN/dE) of all particles produced by DM annihilation.		
    '''
    
    E = logConversion(logE)
    
    int_tot = spec_tot*(E**2)/np.log10(np.e)
    norm = ( conversion(z_dep) )*( trapz(int_tot, logE) )
    #norm = ( conversion(z_dep) )*( 2*mass_dm )

    #int_tot = spec_tot*(E**1)
    #norm = ( conversion(z_dep) )*( trapz(int_tot, E) )
    #norm = ( conversion(z_dep) )*( 2*mass_dm )
                   
    energy_integral = np.zeros( shape=(len(z_inj),len(z_dep)), dtype=np.float64)
    for i in xrange(len(energy_integral)):
        for k in xrange(i,len(energy_integral[i])):
            int_phot = transfer_phot[i,:,k]*spec_phot*(E**2)/np.log10(np.e)
            int_elec = transfer_elec[i,:,k]*spec_elec*(E**2)/np.log10(np.e)
            energy_integral[i][k] = trapz( int_phot + int_elec, logE )

    #energy_integral = np.empty( shape=(len(z_inj),len(z_dep)), dtype=float)
    #for i in xrange(len(energy_integral)):
    #    for k in xrange(i,len(energy_integral[i])):
    #        int_phot = transfer_phot[i,:,k]*spec_phot*(E**1)
    #        int_elec = transfer_elec[i,:,k]*spec_elec*(E**1)
    #        energy_integral[i][k] = trapz( int_phot + int_elec, E )
          
    z_integral = np.zeros_like( z_dep, dtype=np.float64)
    dummy = np.arange(1,len(z_inj)+1)
    for i in xrange(len(z_integral)):
        low = max(i,0)
        #low = i
        integrand = ( conversion(z_inj[low:]) )*energy_integral[i,low:]
        z_integral[i] = trapz( integrand, dummy[low:] )
    
    result = np.empty_like( norm, dtype=np.float64 )
    for i in xrange(len(norm)):
        if norm[i] != 0:
            result[i] = (z_integral[i] / norm[i])
        else:
            result[i] = np.nan
            #result[i] = 0
    
    # Before returning the result-array, we clean the array from NaNs, by interpolating the invalid points
	nanmask = np.isnan(result)
 	#nanmask = np.array([])
    if len(np.flatnonzero(nanmask)) > 0: 
	    #result[nanmask] = np.interp(np.flatnonzero(nanmask), np.flatnonzero(~nanmask), result[~nanmask])	
        pass

    return result

def f_function_decay(logE, z_inj, z_dep, mass_dm, t_dec, transfer_phot, transfer_elec, 
               spec_phot, spec_elec, spec_tot):
    '''
    This is to calculate the f-function given the photon and electron spectrum
    of the particular model, assuming redshift independent DM-annihalation
    i.e. we assume the thermally averaged crossection to be constant. That is
    why we choose alpha=3 as exponent in the conversion factor.   
    
    The inputs are as follows:
    logE:  
        40x1-dim. object containing log10 of the particle kinetic energies 
        (given in units of eV).
    z_inj / z_dep:
        63x1-dim. object containing the log-spaced redshift (z+1) of 
        deposition (or redshift of injection) 
    transfer_phot / transfer_elec:  
        63x40x63-dim. object containing the transfer function for photons 
        (electrons), taking redshift of deposition, Energy at injection and 
        Redshift of injection. 
    spec_phot / spec_elec:  
        Spectrum (dN/dE) of the injected photons as function of Energy 
        (z dependence already factored out).
    spec_tot:
		Total spectrum (dN/dE) of all particles produced by DM annihilation.		
    '''
    
    E = logConversion(logE)
    
    int_tot = spec_tot*(E**2)/np.log10(np.e)
    norm = ( conversion(z_dep,alpha=0) )*( trapz(int_tot, logE) )
    #norm = ( conversion(z_dep,alpha=0) )*( 2*mass_dm )

    #int_tot = spec_tot*(E**1)
    #norm = ( conversion(z_dep,alpha=0) )*( trapz(int_tot, E) )
    #norm = ( conversion(z_dep,alpha=0) )*( 2*mass_dm )
                   
    energy_integral = np.zeros( shape=(len(z_inj),len(z_dep)), dtype=np.float64)
    for i in xrange(len(energy_integral)):
        for k in xrange(i,len(energy_integral[i])):
            int_phot = transfer_phot[i,:,k]*spec_phot*(E**2)/np.log10(np.e)
            int_elec = transfer_elec[i,:,k]*spec_elec*(E**2)/np.log10(np.e)
            energy_integral[i][k] = trapz( int_phot + int_elec, logE )

    #energy_integral = np.empty( shape=(len(z_inj),len(z_dep)), dtype=float)
    #for i in xrange(len(energy_integral)):
    #    for k in xrange(i,len(energy_integral[i])):
    #        int_phot = transfer_phot[i,:,k]*spec_phot*(E**1)
    #        int_elec = transfer_elec[i,:,k]*spec_elec*(E**1)
    #        energy_integral[i][k] = trapz( int_phot + int_elec, E )
          
    z_integral = np.zeros_like( z_dep, dtype=np.float64)
    dummy = np.arange(1,len(z_inj)+1)
    for i in xrange(len(z_integral)):
        low = max(i,0)
        #low = i
        integrand = ( conversion(z_inj[low:],alpha=0) )*np.exp(-time_at_z(z_inj[low:]) / t_dec)*energy_integral[i,low:]
        z_integral[i] = trapz( integrand, dummy[low:] )
    
    result = np.empty_like( norm, dtype=np.float64 )
    for i in xrange(len(norm)):
        if norm[i] != 0:
            result[i] = (z_integral[i] / norm[i])
        else:
            result[i] = np.nan
            #result[i] = 0
    
    # Before returning the result-array, we clean the array from NaNs, by interpolating the invalid points
	nanmask = np.isnan(result)
 	#nanmask = np.array([])
    if len(np.flatnonzero(nanmask)) > 0: 
	    #result[nanmask] = np.interp(np.flatnonzero(nanmask), np.flatnonzero(~nanmask), result[~nanmask])	
        pass

    return result

def log_fit(points,func,xgrid,exponent=1):
    '''
    This is to initialize and call the log(Linear)Interpolator at the same time.
    It is used if the interpolated function needs to be read imediately rather tha stroe it for later use.
    '''
    from interpolator import logInterpolator
    tmp_interpolator = logInterpolator(points, func, exponent)
    #tmp_interpolator = logLinearInterpolator(points, func, exponent)
    out = tmp_interpolator(xgrid)
    return out

def sample_spectrum(input_spec_el, input_spec_ph, input_spec_oth, input_log10E, mass, sampling_log10E, scale='GeV', spec_type='dN/dE'):
    '''
    This method is to sample the input spectra at the energies at which the transfer functions are defined
    '''

    failed = False
    try:
        assert spec_type in ['dN/dE','E.dN/dE']
    except AssertionError:
        print('Unknown type of input spectrum >> %s <<. Valid options are: "dN/dE" and "E.dN/dE"' % spec_type)
        failed = True

    try:
        assert scale in ['ev','keV','MeV','GeV']
    except AssertionError:
        print('Unknown scale of your input >> %s <<. Valid options are: "eV", "keV", "MeV", and "GeV"' % scale )
        failed = True
        scale = 'eV' # Set an abitrary, but valid scale. Since the method is labeled as failed it will return zero anywway

    scale_dict = {
        'GeV': [1e9, 1e-9, 9.],
        'MeV': [1e3, 1e-3, 6.],
        'keV': [1e3, 1e-3, 3.],
        'eV': [1e0, 1e-0, 0.],
    }
    m = mass*scale_dict[scale][0]
    input_log10E += scale_dict[scale][2]*np.ones_like(input_log10E).astype(np.float64)
    # Check if the spectrum is normalized to "integral(E*dN/dE dE) = 2*mass"
    #    If not and the integral is non_zero, rescale the spectrum, if the spectrum is zero (e.g. being in a kinematically forbidden regioern of param space)
    #    return a zero_spectrum
    if spec_type == 'dN/dE':
        factor1 = logConversion(input_log10E)
        factor2 = np.ones_like(input_log10E).astype(np.float64)
    else:
        factor1 = np.ones_like(input_log10E).astype(np.float64)
        factor2 = 1. / logConversion(input_log10E)
    total_dep_energy = trapz( factor1*(input_spec_el+input_spec_ph+input_spec_oth)*scale_dict[scale][1], logConversion(input_log10E) )

    non_zero_spec = ( total_dep_energy > 0.)
    if non_zero_spec and not failed:
        rescaling = total_dep_energy / (2*m)

        out_el = log_fit(input_log10E, factor2*scale_dict[scale][1]*input_spec_el/rescaling, sampling_log10E)
        out_ph = log_fit(input_log10E, factor2*scale_dict[scale][1]*input_spec_ph/rescaling, sampling_log10E)
        out_oth = log_fit(input_log10E, factor2*scale_dict[scale][1]*input_spec_oth/rescaling, sampling_log10E)

        # In the unphysical region E > m the spectrum needs to vanish
        unphysical_region_mask =  (logConversion(sampling_log10E) > m)
        out_el[unphysical_region_mask] = 0.
        out_ph[unphysical_region_mask] = 0.
        out_oth[unphysical_region_mask] = 0.
    else:
        out_el = np.zeros_like(sampling_log10E).astype(np.float64)
        out_ph = np.zeros_like(sampling_log10E).astype(np.float64)
        out_oth = np.zeros_like(sampling_log10E).astype(np.float64)

    return np.array([out_el, out_ph, out_oth])

def finalize(redshift, f_heat, f_lya, f_ionH, f_ionHe, f_lowE):
    '''
    This is to produce the STD-output wich is read by CLASS.
    The structure is as follows.
	 --- The first valid line (i.e. no "#", "%",... or blank) contains the number of lines that follows
     --- Each of the (n) following valid lines contains the table:

     ---->   z_dep | f(z_dep)_heat | f(z_dep)_LyA | f(z_dep)_ionH  | f(z_dep)_ionHe | f(z_dep)_LowE  <--------  
    '''
    first = 2 
    last = len(redshift) - 7
    min_z = 0.
    max_z = 1e4
    print(50*'#')
    print('### This is the standardized output to be read by CLASS.\n### For the correct usage ensure that all other\n### "print(...)"-commands in your script are silenced.')
    print(50*'#'+'\n')
    print('#z_dep\tf_heat\tf_lya\tf_ionH\tf_ionHe\tf_lowE\n\n{:d}\n'.format( (last-first) + 2))
    print('{:.2e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}'.format(min_z,f_heat[first],f_lya[first],f_ionH[first],f_ionHe[first],f_lowE[first]))
    for idx in range(first,last):
	    print('{:.2e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}'.format(redshift[idx],f_heat[idx],f_lya[idx],f_ionH[idx],f_ionHe[idx],f_lowE[idx]))
    print('{:.2e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}\t{:.4e}'.format(max_z,f_heat[last-1],f_lya[last-1],f_ionH[last-1],f_ionHe[last-1],f_lowE[last-1]))
    

