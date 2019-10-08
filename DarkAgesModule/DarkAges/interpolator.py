u"""
.. module:: interpolator
   :synopsis: Definition of the logInterpolator-class(es)
.. moduleauthor:: Patrick Stoecker <stoecker@physik.rwth-aachen.de>

Contains the definition of the classes

- :class:`logInterpolator`: Class to perform a interpolation of a function
  in logarithmic space for a smooth result on a one-dimensional sampling-array.
  (Tries to perform a :code:`cubic` interpolation, if not successful it
  uses a :code:`linear`-interpolation scheme as fallback).
- :class:`logLinearInterpolator`: Similar to the :class:`logInterpolator`-class
  but using the :code:`linear`-interpolation scheme.
- :class:`NDlogInterpolator`: N-dimensional container of instances of the
  :class:`logInterpolator`-class. The underlying points and function to interpolate
  as well as the points to sample at are given as 1D-arrays.

  .. note:: This is not a proper interpolator in N dimensions but a container
	 of 1D-interpolators at a given underlying grid, where the sampling
	 at each of this grid-points needs to be performed for the same array.

"""

from __future__ import absolute_import, division, print_function
from builtins import map, range, object

import numpy as np
from scipy.interpolate import interp1d
from .common import nan_clean
from .__init__ import DarkAgesError

class logInterpolator(object):
	u"""This is a custom class for the purpose of smooth interpolations, if the
	function to interpolate covers many orders of magnitude.

	To do so the interpolation is not done on the function itself, but on its
	logarithm to base 10.

	To improve the quality of the interpolation the function can be skewed
	optionally by multiplying the function by appropriate powers of the abzissa.

	The interplolation is done with the :meth:`interp1d <scipy.interpolate.interp1d>`
	-method of the :class:`scipy.interpolate`-package. For a smooth result
	after backtransformation, the interpolation is done in :code:`cubic`-order.
	The :code:`linear`-order is only used as fallback if the array to interpolate is
	to short for the :code:`cubic` method.

	If the function is only defined at a single point (:code:`len(x) == 1`
	and/or :code:`len(y) == 1`) the interpolator returns the function when
	the point to sample equals the point at which the function is given.
	"""

	def __init__(self, x, y, exponent=1, scale='lin-log'):
		u"""
		Parameters
		----------
		x : :obj:`array-like`
			Array (:code:`shape = (k)`) of points, at which the function
			is defined. **In increasing order**
		y : :obj:`array_like`
		 	Array (:code:`shape = (k)`) with the values of the function
			to interpolate at the points given in :code:`x`
		exponent : :obj:`int`, :obj:`float`, *optional*
			Exponent to specify the powers of :code:`x` mulitplied to
			:code:`y` before the function is transformed into logspace
			(:code:`y -> log10( (x**exponent)*y )`).
			If not specified the default value of :code:`exponent = 1` is taken.
		"""

		if scale not in ['lin-lin','log-lin','lin-log','log-log']:
			raise DarkAgesError('The scale "{0}" is not known. Please choose from ["lin-lin","log-lin","lin-log","log-log".]'.format(scale))
		else:
			self._xscale, self._yscale = scale.split('-')

		self.exponent = exponent
		self._lower = float(x[0])
		self._upper = float(x[-1])
		if len(x) != len(y):
			raise DarkAgesError('Array of x-values does not correspond to the array of y-values (different sizes)')
		if len(x) < 2:
			self._only_one_point = True
			self._single_x_value = float(x[0])
			self._single_y_value = float(y[0])
			return None
		else:
			self._only_one_point = False

		func_copy = np.zeros_like(y, dtype=np.float64)
		for idx in range(len(y)):
			if (y[idx] < 0) or (y[idx] != y[idx]):
				func_copy[idx] = np.nan
				#func_copy[idx] = 0
			else:
				func_copy[idx] = y[idx]

		nan_mask = func_copy == func_copy

		if self._yscale == 'log':
			fitfunc = np.log1p( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
		else:
			fitfunc = func_copy[nan_mask] * (x[nan_mask]**self.exponent)

		if self._xscale == 'log':
			points = np.log1p(x[nan_mask])
		else:
			points = x[nan_mask]

		self._valid = True
		if len(func_copy[nan_mask]) > 4:
			self._fit = interp1d(points, fitfunc, bounds_error=False, fill_value=np.nan, kind='cubic')
		elif len(func_copy[nan_mask]) >= 2:
			self._fit = interp1d(points, fitfunc, bounds_error=False, fill_value=np.nan, kind='linear')
		else:
			self._valid = False
		return None

	def __call__(self, xgrid):
		u"""
		Parameters
		----------
		xgrid : :obj:`array-like`
			Array (:code:`shape = (l)`) with the points at which the function
			should be interpolated.

		Returns
		-------
		:obj:`array-like`
			Array (:code:`shape = (l)`) with the interpolated values of the
			original function (:code:`y`). For points in :code:`xgrid`
			outside the bounds of :code:`x`, :math:`0` is returned.
		"""

		def dummy(single_point):
			if self._only_one_point:
				if abs((single_point - self._single_x_value)/single_point) <= 1e-5:
					return self._single_y_value
				else:
					return 0.
			else:
				if self._valid:
					if self._xscale == 'log':
						xval = np.log1p(single_point)
					else:
						xval = single_point
					if self._yscale == 'log':
						return np.expm1(self._fit(xval)) / (single_point**self.exponent)
					else:
						return (self._fit(xval)) / (single_point**self.exponent)
				else:
					return 0.

		out = nan_clean( np.vectorize(dummy).__call__(xgrid) )
		return out

	def get_upper(self):
		u"""Returns the upper bound of the domain of the function.

		Parameters
		----------
		None

		Returns
		-------
		:obj:`float`
			Upper bound of the domain (:code:`x[-1]`)
		"""
		return self._upper

	def get_lower(self):
		u"""Returns the lower bound of the domain of the function.

		Parameters
		----------
		None

		Returns
		-------
		:obj:`float`
			Lower bound of the domain (:code:`x[0]`)
		"""

		return self._lower

class NDlogInterpolator(object):
	u"""This is a custom class with the purpose to act as a container for
	several instances of the
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`-class
	of functions defined on the same points :code:`x` and interpolated at the
	same points :code:`xgrid` afterwards.

	This is for example used of you have the spectra of electrons and positrons,
	photons, and other particles which are given as a table at the
	same energies and you want to sample each of this three spectra at the same
	points.
	Instead of looping over the three instances of the
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`-class
	to and calling
	:meth:`logInterpolator.__call__ <DarkAges.interpolator.logInterpolator.__call__>`
	for each of the spectra, a single-line call to
	:meth:`NDlogInterpolator.__call__ <DarkAges.interpolator.NDlogInterpolator.__call__>`
	does this loop automatically.

	.. note:: There is no gain in runtime with using this class. The purpose
	   is more to provide a shortcut for a cleaner code.

	For further information on the interpolation, see
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`
	"""

	def __init__(self, x, array_of_y, exponent=1, scale='lin-log'):
		u"""
		Parameters
		----------
		x : :obj:`array-like`
			Array (:code:`shape = (k)`) of points, at which the function
			is defined. **In increasing order**
		array_of_y : :obj:`array_like`
		 	Array (:code:`shape = (k, m1,...,mn)`) with the values of the functions
			to interpolate at the points given in :code:`x`. First dimension refers
			to the points at which the functions are defined. All other dimensions
			refer to the shape of the multidimensional grid of the functions to
			interpolate
		exponent : :obj:`int`, :obj:`float`, *optional*
			Exponent to specify the powers of :code:`x` mulitplied to
			:code:`array_of_y` (its first dimension) before the function is
			transformed into logspace (:code:`array_of_y -> log10( (x**exponent)*array_of_y )`).
			If not specified the default value of :code:`exponent = 1` is taken.
		"""
		if scale not in ['lin-lin','log-lin','lin-log','log-log']:
			raise DarkAgesError('The scale "{0}" is not known. Please choose from ["lin-lin","log-lin","lin-log","log-log".]'.format(scale))
		self._ndim = array_of_y.ndim
		self._shape = array_of_y.shape[1:]
		self._lower = float(x[0])
		self._upper = float(x[-1])
		if self._ndim == 1:
			self.interpolator = logInterpolator(x, array_of_y, exponent=exponent)
			#self.interpolator = logLinearInterpolator(x, array_of_y, exponent=exponent)
		elif self._ndim > 1:
			self.interpolator = np.empty(shape=self._shape, dtype=logInterpolator)
			rolled_y = np.rollaxis(array_of_y,0,self._ndim)
			for multi_idx in np.ndindex(self.interpolator.shape):
				#self.interpolator[multi_idx] = logInterpolator(x[:], rolled_y[multi_idx], exponent=exponent, , scale=scale)
				self.interpolator[multi_idx] = logLinearInterpolator(x[:], rolled_y[multi_idx], exponent=exponent, scale=scale)

	def __call__(self, xgrid):
		u"""
		Parameters
		----------
		xgrid : :obj:`array-like`
			Array (:code:`shape = (l)`) with the points at which the functions
			should be interpolated.

		Returns
		-------
		:obj:`array-like`
			Array (:code:`shape = (l,m1,...,mn)`) with the interpolated values of the
			original functions (:code:`array_of_y`). For points in :code:`xgrid`
			outside the bounds of :code:`x`, :math:`0` is returned.
		"""
		def dummy(single_point):
			if self._ndim == 1:
				return self.interpolator.__call__(single_point)
			elif self._ndim > 1:
				ret = np.zeros(shape=self._shape, dtype=np.float64)
				for multi_idx in np.ndindex(self.interpolator.shape):
					ret[multi_idx] = nan_clean( self.interpolator[multi_idx].__call__(single_point) )
				return ret

		out = np.vectorize(dummy, otypes=[np.ndarray]).__call__(xgrid)
		try:
			out = np.array(list(map(list, out)), dtype=np.float64)
		except TypeError:
			out = np.array(out, dtype=np.float64)
		return out

	def get_upper(self):
		u"""Returns the upper bound of the domain of the function.
		See :meth:`logInterpolator.get_upper <DarkAges.interpolator.logInterpolator.get_upper>`
		"""
		return self._upper

	def get_lower(self):
		u"""Returns the lower bound of the domain of the function.
		See :meth:`logInterpolator.get_lower <DarkAges.interpolator.logInterpolator.get_lower>`
		"""
		return self._lower

class logLinearInterpolator(object):
	u"""This is a custom class for the purpose of smooth interpolations, if the
	function to interpolate covers many orders of magnitude.

	It is defined similar to the
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`-class with the only
	difference that the interpolation is done with the :code:`linear`-method.

	It has the same attributes and methods as the
	:class:`logInterpolator <DarkAges.interpolator.logInterpolator>`-class. Please
	look there for the definition of the parameters and returned data types.
	"""

	def __init__(self, x, y, exponent=1, scale='lin-log'):
		if scale not in ['lin-lin','log-lin','lin-log','log-log']:
			raise DarkAgesError('The scale "{0}" is not known. Please choose from ["lin-lin","log-lin","lin-log","log-log".]'.format(scale))
		else:
			self._xscale, self._yscale = scale.split('-')

		self.exponent = exponent
		self._lower = float(x[0])
		self._upper = float(x[-1])
		if len(x) != len(y):
			raise DarkAgesError('Array of x-values does not correspond to the array of y-values (different sizes)')
		if len(x) < 2:
			self._only_one_point = True
			self._single_x_value = float(x[0])
			self._single_y_value = float(y[0])
			return None
		else:
			self._only_one_point = False

		func_copy = np.zeros_like(y, dtype=np.float64)
		for idx in range(len(y)):
			if (y[idx] <= 0) or (y[idx] != y[idx]):
				func_copy[idx] = np.nan
				#func_copy[idx] = 0
			else:
				func_copy[idx] = y[idx]

		nan_mask = func_copy == func_copy

		self._valid = True
		if len(func_copy[nan_mask]) >= 2:
			if self._yscale == 'log':
				fitfunc = np.log1p( func_copy[nan_mask] * (x[nan_mask]**self.exponent) )
			else:
				fitfunc = func_copy[nan_mask] * (x[nan_mask]**self.exponent)
			if self._xscale == 'log':
				points = np.log1p(x[nan_mask])
			else:
				points = x[nan_mask]
			self._fit = interp1d(points, fitfunc, bounds_error=False, fill_value=np.nan, kind='linear')
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
					if self._xscale == 'log':
						xval = np.log1p(single_point)
					else:
						xval = single_point
					if self._yscale == 'log':
						return np.expm1(self._fit(xval)) / (single_point**self.exponent)
					else:
						return (self._fit(xval)) / (single_point**self.exponent)
				else:
					return 0.

		out = nan_clean( np.vectorize(dummy).__call__(xgrid) )
		return out

	def get_upper(self):
		return self._upper

	def get_lower(self):
		return self._lower
