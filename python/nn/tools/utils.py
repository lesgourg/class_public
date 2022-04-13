import functools
from itertools import islice
from contextlib import contextmanager
import time
import torch

# some constants

# speed of light
C = 2997.92458

def chunk(it, size):
    it = iter(it)
    return iter(lambda: tuple(islice(it, size)), ())

class dotdict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

@contextmanager
def timing(message=None, fmt="{:.3f}"):
    if message is not None:
        print(message)
    start = time.time()
    yield
    end = time.time()
    duration = end - start
    print("done; took", fmt.format(duration), "seconds.")

def transpose_list_of_dicts(lst):
    """
    Given a list of dicts with shared keys, return a single dict
    whose values are lists of the values.

    Inverse of `transpose_dict_of_lists`.
    """
    keys = [frozenset(dct) for dct in lst]
    all_keys_equal = len(set(keys)) <= 1
    assert all_keys_equal
    keys = keys[0]

    return {key: [dct[key] for dct in lst] for key in keys}


def transpose_dict_of_lists(dct):
    """
    Given a  dict of lists, return a list of dicts.

    Inverse of `transpose_list_of_dicts`.
    """
    size = len(next(iter(dct.values())))
    #list of parameters which only appear once
    keys_once = [key for key in dct.keys() if len(dct[key])==1]
    keys_size = [key for key in dct.keys() if len(dct[key])==size]

    out_dict = []

    for i in range(size):
        _ = {key: dct[key][i] for key in keys_size}
        for key in keys_once:
            _[key] = dct[key][0]
        out_dict.append(_)

    return out_dict


def powerspace(start, stop, k, *args, **kwargs):
    """Returns an array such that arr**(1/k) is equivalent to linspace"""
    start = start ** (1/k)
    stop = stop ** (1/k)

    import torch
    return torch.linspace(start, stop, *args, **kwargs)**k

def torch_gradient(y, x):
    """ Compute the gradient of y w.r.t. x along the last axis, i.e. y has shape (dim0, dim1, ..., len(x)) """
    result = torch.empty(y.shape, device=y.device)
    result[..., 1:-1] = (y[..., 2:] - y[..., :-2]) / (x[..., 2:] - x[..., :-2])
    # extrapolate linearly at the boundaries (not optimal, but whatever...)
    result[..., 0] = result[..., 1] - (result[..., 2] - result[..., 1]) / (x[..., 2] - x[..., 1]) * (x[..., 1] - x[..., 0])
    result[..., -1] = result[..., -2] + (result[..., -2] - result[..., -3]) / (x[..., -2] - x[..., -3]) * (x[..., -2] - x[..., -3])
    return result


def extrapolate_left(k, S, k_):
    """
    Given an array `k` and a 2D function `S` of shape (n, len(k)),
    return a tuple of a new k array and S array,
    in which S has been linearly extrapolated to those values
    out of `k_` which are below the smallest value in `k`.
    """
    import numpy as np
    k_ = k_[k_ < k[0]]
    slope = (S[:, 1] - S[:, 0]) / (k[1] - k[0])
    dk = k_ - k[0]
    S_ = S[:, [0]] + slope[:, None] * dk[None, :]
    k_ = np.concatenate((k_, k))
    S_ = np.concatenate((S_, S), axis=1)
    return k_, S_

def extrapolate_right(k, S, k_):
    """
    Like `extrapolate_left`, but on the right side
    """
    import numpy as np
    k_ = k_[k_ > k[-1]]
    slope = (S[:, -1] - S[:, -2]) / (k[-1] - k[-2])
    dk = k_ - k[-1]
    S_ = S[:, [-1]] + slope[:, None] * dk[None, :]
    k_ = np.concatenate((k, k_))
    S_ = np.concatenate((S, S_), axis=1)
    return k_, S_

def extrapolate_left_pow(k, S, k_):
    """
    Given an array `k` and a 2D function `S` of shape (n, len(k)),
    return a tuple of a new k array and S array,
    in which S has been extrapolated using a power law to those values
    out of `k_` which are below the smallest value in `k`.
    """
    import numpy as np
    k_ = k_[k_ < k[0]]
    b = np.log(S[:, 1] / S[:, 0]) / np.log(k[1] / k[0])
    a = S[:, 0] / k[0]**b
    S_ = a[:, None] * k_[None, :]**b[:, None]
    k_ = np.concatenate((k_, k))
    S_ = np.concatenate((S_, S), axis=1)
    return k_, S_

def extrapolate_right_pow(k, S, k_):
    """
    Like `extrapolate_left_pow`, but on the right side
    """
    import numpy as np
    k_ = k_[k_ > k[-1]]
    b = np.log(S[:, -1] / S[:, -2]) / np.log(k[-1] / k[-2])
    a = S[:, -1] / k[-1]**b
    S_ = a[:, None] * k_[None, :]**b[:, None]
    k_ = np.concatenate((k, k_))
    S_ = np.concatenate((S, S_), axis=1)
    return k_, S_

def extrapolate(k, S, k_):
    k, S = extrapolate_left(k, S, k_)
    assert len(k) == S.shape[1]
    k, S = extrapolate_right(k, S, k_)
    assert len(k) == S.shape[1]
    return k, S

def extrapolate_pow(k, S, k_):
    k, S = extrapolate_left_pow(k, S, k_)
    assert len(k) == S.shape[1]
    k, S = extrapolate_right_pow(k, S, k_)
    assert len(k) == S.shape[1]
    return k, S

@contextmanager
def timeit(label, sync=False):
    start = time.perf_counter()
    yield
    if sync:
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - start
    print("== TIMING == Action `{}` took {:.3e}s.".format(label, elapsed))

class Timer:
    """
    Simple help for performance measurements.
    """
    def __init__(self):
        self._start = {}
        self._end = {}
        self._times = {}

    def start(self, name):
        if name in self._start:
            print("WARNING: Overwriting measurement {}".format(name))
        self._start[name] = time.perf_counter()

    def end(self, name):
        if name not in self._start:
            raise ValueError(
               "Measurement '{}' has not started; cannot end!".format(name)
               )
        if name in self._end:
            print("WARNING: Overwriting measurement {}".format(name))
        self._end[name] = time.perf_counter()
        self._times[name] = self._end[name] - self._start[name]

    @property
    def times(self):
        return self._times

    def __getitem__(self, name):
        return self._times[name]

    def __setitem__(self, name, value):
        self._times[name] = value

    @contextmanager
    def timeit(self, name):
        self.start(name)
        yield
        self.end(name)

