'''
This module provides the customized (1D)-interpolator classes.
	- logInterpolator (using 'cubic' if more than 2 data points, 'linear' if only two)
	- logLinearInterpolator (using 'linear' for more or equal than two data points)

For a smooth interpolation the points are transformed into log-space 
(optionally the points are multiplied by powers of the abzissa beforehand) 
afterwards transformed backwards
'''

import numpy as np
from scipy.interpolate import interp1d
import common as misc

class logInterpolator(object):
    def __init__(self, x, y, exponent=1):
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
            out = misc.nan_clean( tmp_out )
        return out

class logLinearInterpolator(object):
    def __init__(self, x, y, exponent=1):
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
            out = misc.nan_clean( tmp_out )
        return out

