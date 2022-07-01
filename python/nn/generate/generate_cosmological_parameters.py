import numpy as np
from pyDOE import lhs

def sample_cosmological_parameters(domain, count):
    """
    Sample cosmological parameters from `domain` (which is specified as a dict
    names -> (min, max)) using Latin Hypercube Sampling.
    The parameters are returned as a dict mapping the keys from `domain` to numpy arrays.
    """
    # array of shape (#samples, #dimensions) with values from (0, 1)
    sampling = lhs(n=len(domain), samples=count, criterion="center")
    return {key: lo + (hi - lo) * column for (key, (lo, hi)), column in zip(domain.items(), sampling.T)}

def write_to_h5(samples, path):
    # Store data
    import h5py as h5
    with h5.File(output_filepath, "w") as f:
        for quantity, array in samples.items():
            f.create_dataset(quantity, dtype="f", data=array)
