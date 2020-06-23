import numpy as np

from .data_providers.providers import CLASSDataProvider

class PredictorCache:
    def __init__(self, raw_inputs, transformed_inputs):
        self.raw_inputs = raw_inputs
        self.transformed_inputs = transformed_inputs

    def get_raw_inputs(self, names, tau_mask=None):
        return self._get_inputs(self.raw_inputs, names, tau_mask=tau_mask)

    def get_transformed_inputs(self, names, tau_mask=None):
        return self._get_inputs(self.transformed_inputs, names, tau_mask=tau_mask)

    def _get_inputs(self, source, names, tau_mask=None):
        # TODO names
        def cut(arr):
            if tau_mask is None:
                return arr
            else:
                return arr[tau_mask, ...]

        return {k: cut(source[k]) for k in names}




