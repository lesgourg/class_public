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
from .common import nan_clean, print_warning

class logInterpolator(object):
	def __init__(self, x, y, exponent=1):
		self._valid = True
		self.exponent = exponent
		self._lower = x[0]
		self._upper = x[-1]
		try:
			assert len(x) == len(y)
		except AssertionError:
			print_warning('Array of x-values does not correspond to the array of y-values (different sizes)')
			self._valid = False
			return None
		if len(x) < 2:
			self._only_one_point = True
			self._single_x_value = x
			self._single_y_value = y
			return None
		else:
			self._only_one_point = False

		func_copy = np.zeros_like(y, dtype=np.float64)
		for idx in xrange(len(y)):
			if (y[idx] <= 0) or (y[idx] != y[idx]):
				func_copy[idx] = np.nan
				#func_copy[idx] = 0
			else:
				func_copy[idx] = y[idx]

		nan_mask = func_copy == func_copy

		if len(func_copy[nan_mask]) > 4:
			log_f = np.log( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
			self.log_f_fit = interp1d(x[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='cubic')
		elif len(func_copy[nan_mask]) >= 2:
			log_f = np.log( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
			self.log_f_fit = interp1d(x[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='linear')
		else:
			self._valid = False
		return None

	def __call__(self, xgrid):
		def dummy(single_point):
			if self._only_one_point:
				if abs((single_point - self._single_x_value)/single_point) <= 1e-5:
					return self._single_y_value
				else:
					return 0.
			else:
				if self._valid:
					return np.e**(self.log_f_fit(single_point)) / (single_point**self.exponent)
				else:
					return 0.

		out = nan_clean( np.vectorize(dummy).__call__(xgrid) )
		return out

	def get_upper(self):
		return self._upper

	def get_lower(self):
		return self._lower

class NDlogInterpolator(object):
	def __init__(self, x, array_of_y, exponent=1):
		self._ndim = array_of_y.ndim
		self._shape = array_of_y.shape[1:]
		self._lower = x[0]
		self._upper = x[-1]
		if self._ndim == 1:
			self.interpolator = logInterpolator(x, array_of_y, exponent=exponent)
			#self.interpolator = logLinearInterpolator(x, array_of_y, exponent=exponent)
		elif self._ndim > 1:
			self.interpolator = np.empty(shape=self._shape, dtype=logInterpolator)
			rolled_y = np.rollaxis(array_of_y,0,self._ndim)
			for multi_idx in np.ndindex(self.interpolator.shape):
				#self.interpolator[multi_idx] = logInterpolator(x[:], rolled_y[multi_idx], exponent=exponent)
				self.interpolator[multi_idx] = logLinearInterpolator(x[:], rolled_y[multi_idx], exponent=exponent)

	def __call__(self, xgrid):
		def dummy(single_point):
			if self._ndim == 1:
				return self.interpolator.__call__(single_point)
			elif self._ndim > 1:
				ret = np.zeros(shape=self._shape, dtype=np.float64)
				for multi_idx in np.ndindex(self.interpolator.shape):
					ret[multi_idx] = nan_clean( self.interpolator[multi_idx].__call__(single_point) )
				return ret

		out = np.vectorize(dummy, otypes=[np.ndarray]).__call__(xgrid)
		out = np.array(map(list, out), dtype=np.float64)
		return out

	def get_upper(self):
		return self._upper

	def get_lower(self):
		return self._lower

class logLinearInterpolator(object):
	def __init__(self, x, y, exponent=1):
		self._valid = True
		self.exponent = exponent
		self._lower = x[0]
		self._upper = x[-1]
		try:
			assert len(x) == len(y)
		except AssertionError:
			print_warning('Array of x-values does not correspond to the array of y-values (different sizes)')
			self._valid = False
			return None
		if len(x) < 2:
			self._only_one_point = True
			self._single_x_value = x
			self._single_y_value = y
			return None
		else:
			self._only_one_point = False

		func_copy = np.zeros_like(y, dtype=np.float64)
		for idx in xrange(len(y)):
			if (y[idx] <= 0) or (y[idx] != y[idx]):
				func_copy[idx] = np.nan
				#func_copy[idx] = 0
			else:
				func_copy[idx] = y[idx]

		nan_mask = func_copy == func_copy

		if len(func_copy[nan_mask]) >= 2:
			log_f = np.log( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
			self.log_f_fit = interp1d(x[nan_mask], log_f, bounds_error=False, fill_value=np.nan, kind='linear')
		else:
			self._valid = False
		return None

	def __call__(self, xgrid):
		def dummy(single_point):
			if self._only_one_point:
				if abs((single_point - self._single_x_value)/single_point) <= 1e-5:
					return self._single_y_value
				else:
					return 0.
			else:
				if self._valid:
					return np.e**(self.log_f_fit(single_point)) / (single_point**self.exponent)
				else:
					return 0.

		out = nan_clean( np.vectorize(dummy, excluded=['isvalid']).__call__(xgrid) )
		return out

	def get_upper(self):
		return self._upper

	def get_lower(self):
		return self._lower
