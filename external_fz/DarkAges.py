# -*- coding: utf-8 -*-
"""
@author: patrick
"""

from scipy.integrate import trapz
from scipy.interpolate import interp1d, interp2d
import dill
import numpy as np

def logConversion(LogNumber):
    return 10**LogNumber
    
def H(redshift):
    '''
    Redshift dependence of H (assuming only matter domination). 
    The normalization is dropped, since it cancels anyway
    in the calculation of f.
    '''
    return (redshift)**(3./2.)

def conversion( redshift, alpha ):
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

def f_function(logE, z_inj, z_dep, mass_dm, transfer_phot, transfer_elec, 
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
    norm = ( conversion(z_dep,3) )*( trapz(int_tot, logE) )
    #norm = ( conversion(z_dep,3) )*( 2*mass_dm )

    #int_tot = spec_tot*(E**1)
    #norm = ( conversion(z_dep,3) )*( trapz(int_tot, E) )
    #norm = ( conversion(z_dep,3) )*( 2*mass_dm )
                   
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
        integrand = ( conversion(z_inj[low:],3) )*energy_integral[i,low:]
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

def log_fit(points,func,xgrid):
    #try:
    #    assert len(points) == len(func)
    #except AssertionError:
    #    print 'Array of x-values does not correspond to the array of y-values (different sizes)'
    #    return np.zeros_like(xgrid)
    #
    #f_copy = np.zeros_like(func, dtype=np.float64)
    #for idx in xrange(len(func)):
    #    if (func[idx] <= 0) or (func[idx] != func[idx]):
    #        f_copy[idx] = np.nan
    #        #f_copy[idx] = 0
    #    else:
    #        f_copy[idx] = func[idx]
    #
    #nan_mask = f_copy == f_copy
    #
    #if len(f_copy[nan_mask]) > 0:
    #    log_f = np.log( f_copy[nan_mask] * (points[nan_mask]**2) )
    #    log_f_fit = interp1d(points[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='cubic')
    #    out = nan_clean( np.e**(log_f_fit(xgrid)) / (xgrid**2) )
    #else:
    #    out = np.zeros_like(xgrid, dtype=np.float64)
    #
    #return out
    
    tmp_interpolator = logInterpolator(points, func, 1)
    #tmp_interpolator = logLinearInterpolator(points, func, 1)
    out = tmp_interpolator(xgrid)
    return out


class logInterpolator(object):
    def __init__(self, x, y, exponent):
        self.valid = True
        self.exponent = exponent
        try:
            assert len(x) == len(y)
        except AssertionError:
            print 'Array of x-values does not correspond to the array of y-values (different sizes)'
            self.valid = False
            return None
    
        func_copy = np.zeros_like(y, dtype=np.float64)
        for idx in xrange(len(y)):
            if (y[idx] <= 0) or (y[idx] != y[idx]):
                func_copy[idx] = np.nan
                #func_copy[idx] = 0
            else:
                func_copy[idx] = y[idx]
    
        nan_mask = func_copy == func_copy

        #print '%i out of %i values are masked.' % (len(func_copy[~nan_mask]),len(func_copy))

        if len(func_copy[nan_mask]) >= 3:
            log_f = np.log( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
            self.log_f_fit = interp1d(x[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='cubic')
        elif len(func_copy[nan_mask]) >= 2:
            log_f = np.log( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
            self.log_f_fit = interp1d(x[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='linear')
        else:
            self.valid = False
        return None

    def __call__(self, xgrid):
        if xgrid.shape == ():
			xgrid = np.asarray([xgrid])
        if not self.valid:
            out = np.zeros_like(xgrid, dtype=np.float64)
        else:
            tmp_out = np.e**(self.log_f_fit(xgrid)) / (xgrid**self.exponent)
            out = nan_clean( tmp_out )
        return out

class logLinearInterpolator(object):
    def __init__(self, x, y, exponent):
        self.valid = True
        self.exponent = exponent
        try:
            assert len(x) == len(y)
        except AssertionError:
            print 'Array of x-values does not correspond to the array of y-values (different sizes)'
            self.valid = False
            return None
    
        func_copy = np.zeros_like(y, dtype=np.float64)
        for idx in xrange(len(y)):
            if (y[idx] <= 0) or (y[idx] != y[idx]):
                func_copy[idx] = np.nan
                #func_copy[idx] = 0
            else:
                func_copy[idx] = y[idx]
    
        nan_mask = func_copy == func_copy

        #print '%i out of %i values are masked.' % (len(func_copy[~nan_mask]),len(func_copy))

        if len(func_copy[nan_mask]) >= 2:
            log_f = np.log( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
            self.log_f_fit = interp1d(x[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='linear')
        else:
            self.valid = False
        return None

    def __call__(self, xgrid):
        if xgrid.shape == ():
			xgrid = np.asarray([xgrid])
        if not self.valid:
            out = np.zeros_like(xgrid, dtype=np.float64)
        else:
            tmp_out = np.e**(self.log_f_fit(xgrid)) / (xgrid**self.exponent)
            out = nan_clean( tmp_out )
        return out

class logGridInterpolator(object):
    def __init__(self, x, y, z):
        self.valid = True
        try:
            assert len(x)*len(y) == z.shape[0]*z.shape[1]
        except AssertionError:
            print 'Arrays of x-values and y-values do not correspond to the 2D-grid of z-values (different sizes)'
            self.valid = False
            return None

        x_list = np.repeat(x, len(y)).astype(np.float64)
        y_list = np.tile(y, len(x)).astype(np.float64)
        z_list = np.zeros(shape=(len(x)*len(y)), dtype=np.float64)
        try:
            assert len(x_list) == len(y_list) and  len(x_list) == len(y_list)
        except AssertionError:
            print 'Arrays of x-values and y-values do not correspond to the 2D-grid of z-values (different sizes)'
            self.valid = False
            return None
        for idx1 in xrange(len(x)):
            for idx2 in xrange(len(y)):
                tot_idx = len(y)*idx1+idx2
                if (z[idx1,idx2] <= 0) or (z[idx1,idx2] != z[idx1,idx2]):
                    z_list[tot_idx] = np.nan
                else:
                    z_list[tot_idx] = z[idx1,idx2]
    
        nan_mask = z_list == z_list

        log_z = np.log( z_list[nan_mask] )
        self.log_z_fit = interp2d(x_list[nan_mask],y_list[nan_mask], log_z, bounds_error=False, fill_value=np.nan, kind='cubic')
        return None

    def __call__(self, xgrid, ygrid):
        if xgrid.shape == ():
            xgrid = np.asarray([xgrid])
        if ygrid.shape == ():
            ygrid = np.asarray([ygrid])
        out = np.zeros(shape=(len(xgrid),len(ygrid)), dtype=np.float64)
        if self.valid:
            for idx1 in xrange(out.shape[0]):
                for idx2 in xrange(out.shape[1]):
                    tmp_out = np.e**(self.log_z_fit(xgrid[idx1],ygrid[idx2]))
                    out[idx1,idx2] = nan_clean( tmp_out )
        return out


class transfer(object):
    def __init__(self, infile):
        print 'Initializing the transfer functions'
        data = np.genfromtxt(infile, unpack=True, usecols=(0,1,2,3,4), dtype=np.float64 )
        self.z_injected = np.unique(data[2]).astype(np.float64)
        self.z_deposited = np.unique(data[0]).astype(np.float64)
        self.log10E = np.unique(data[1]).astype(np.float64)
        l1 = len(self.z_deposited)
        l2 = len(self.log10E)
        l3 = len(self.z_injected)
        self.transfer_phot = data[4].reshape(l1,l2,l3).astype(np.float64)
        self.transfer_elec = data[3].reshape(l1,l2,l3).astype(np.float64)   
        print 'Initializing the transfer functions: Done!'
 
class model(object):
    def __init__(self,el_spec,ph_spec,oth_spec,m):
        self.mass = m
        self.photon_spec = ph_spec
        self.electron_spec = el_spec
        self.total_spec = ph_spec+el_spec+oth_spec
    
    def calc_f(self, transfer_instance):
        if not isinstance(transfer_instance, transfer):
            print 'Warning!!! You did not include a proper instance of the class "transfer"'
            return -1
        else:
            red = transfer_instance.z_deposited
            f_func = f_function(transfer_instance.log10E, transfer_instance.z_injected, 
                                transfer_instance.z_deposited, self.mass,
                                transfer_instance.transfer_phot, 
                                transfer_instance.transfer_elec, 
                                self.photon_spec, self.electron_spec, self.total_spec)
            return np.array([red, f_func], dtype=np.float64)
    
    def save_f(self,transfer_instance, filename):
        f_function = self.calc_f(transfer_instance)
        file_out = open(filename, 'w')
        file_out.write('#z_dep\tf(z)')
        for i in range(len(f_function[0])):
            file_out.write('\n%.4g\t%.8g'%(f_function[0,i],f_function[1,i]))
        file_out.close()
        print 'Saved effective f(z)-curve under "%s"'%filename
