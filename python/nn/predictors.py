import multiprocessing
import functools
import time

from .data_providers import CLASSDataProvider
# from plotting.plot_source_function import plot_source_function
# import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import os
from time import perf_counter
import scipy.interpolate
import math
import torch
import matplotlib.pyplot as plt
from .models import ALL_NETWORK_CLASSES
from .data_providers import current_transformer

def lowest_k_index(k, k_min):
    idx = np.searchsorted(k, k_min)
    return max(idx - 1, 0)

class BasePredictor:
    def __init__(self, cosmo, input_transformer, target_transformer, k):
        self.cosmo = cosmo
        self.k = k
        self.tau_size = 0
        self.k_size = 0

        self.input_transformer = input_transformer
        self.target_transformer = target_transformer

        self.k_min_class = self.cosmo.k_min()
        self.k_min_idx = lowest_k_index(self.k, self.k_min_class)

        self.reset_cache()
        self.reset_times()

    def cleanup(self):
        # Called by classy in Class.struct_cleanup()
        self.reset_cache()
        self.reset_times()

    def reset_cache(self):
        self.cache = {}

    def reset_times(self):
        # NOTE: use explicit initializiation instead of defaultdict
        # to catch bad access on non existing keys
        self.times = {
            # time spent in `predict()` and children
            "predictor.predict":                    0.0,
            # time spent only in network evaluation
            "neural network evaluation":            0.0,
            "neural network input transformation":  0.0,
            "neural network output transformation": 0.0,
            "save source function in cache": 0.0,
            "interpolate k_min": 0.0,
            "build DataProvider": 0.0,
            "single predictions complete":0.0,
            "model load": 0.0,
            "cache load": 0.0,
            "rule application": 0.0,
            "float to double transformation": 0.0,
        }
        self.time_prediction_per_network = {}

    def get_inputs(self, cosmo, tau, selection, mask=None, cached=True):
        """
        This method will be called by `predict_many` to obtain the necessary
        inputs for the NN evaluation from the instance `cosmo` of CLASS
        for the given sampling of `tau`.
        This will return a dictionary of 'raw inputs' (i.e. having physically meaningful values)
        and a dictionary of 'transformed inputs' (i.e. the normalized version used for NN evaluation).
        """
        raise DeprecationWarning("get_inputs is deprecated!")

        # if self.raw_inputs and cached:
        #     assert self.transformed_inputs
        #     all_cached = all(item in self.raw_inputs for item in selection)
        #     all_cached = all_cached and all(item in self.transformed_inputs for item in selection)
        #     if all_cached:
        #         # make sure that we know the inputs at all requested tau values
        #         # by checking that all tau values are in self.raw_inputs["tau"]
        #         all_tau = np.sum(self.raw_inputs["tau"] == tau) == len(tau)
        #         assert all_tau
        #         tau_subset = self.raw_inputs["tau"] == tau
        #         raw_inputs_subset = {key: value[tau_subset] for key, value in self.raw_inputs}
        #         tf_inputs_subset = {key: value[tau_subset] for key, value in self.tf_inputs}
        #         # return self.raw_inputs, self.transformed_inputs
        #         return raw_inputs_subset, tf_inputs_subset

        self.create_provider_if_not_exists(cosmo)
        raw_inputs = self.provider.get_inputs(k=self.k, tau=tau, input_selection=selection)
        assert all(item in raw_inputs for item in selection)
        transformed_inputs = self.input_transformer.transform_inputs(raw_inputs)
        assert all(item in transformed_inputs for item in selection)

        return raw_inputs, transformed_inputs

    def get_limit(self, quantity, raw_inputs):
        """
        This method returns the k -> 0 limit of the source function named `quantity`.
        """
        if quantity in ("t0", "t0_reco", "t0_sw", "t0_reco_no_isw", "t0_reio", "t0_reio_no_isw"):
            return raw_inputs["g"] / 3
        elif quantity == "phi_plus_psi":
            return 4. / 3.
        elif quantity == "delta_m":
            return -6. / 5.
        elif quantity == "delta_cb":
            return -6. / 5.
        else:
            return 0

    def predict(self, source_array, quantity, tau, provider=None, cache=None):
        """
        Get the network prediction for the given `quantity` and `tau` array.
        The result will be sampled on the `k` array obtained from `self.get_k()` and
        thus have the shape `(len(self.get_k()), len(tau))`
        """

        source_array = self._predict(quantity, tau, provider, cache, source_array)

        return source_array

    def _predict(self, quantity, tau, source_array):
        """
        Predict source function for given quantity.
        Will be implemented by child classes.
        """
        raise NotImplementedError

    def _all_network_input_names(self):
        """
        Return list of all network inputs required.
        Will be implemented by child classes.
        """
        raise NotImplementedError

    def get_k(self):
        k = self.k

        # NOTE: THIS IS VERY IMPORTANT! WE MAY NOT PASS K MODES < K_MIN_CLASS TO CLASS
        # OTHERWISE BAD THINGS WILL HAPPEN!!!
        k_result = np.concatenate((
            np.array([self.k_min_class]),
            k[k > self.k_min_class]
        ))
        return k_result

    def predict_many(self, quantities, tau, sources_array = None):
        """
        Predict the source functions whose names are given as the list `quantities`.
        This will return a numpy array of shape (len(quantities), len(k) + 1, len(tau)).
        The 2nd size is len(k) + 1 instead of len(k) because this function adds another
        row to the S array corresponding to k -> 0.
        This is needed so that CLASS does not have to extrapolate for low k.
        """

        # load tau size
        self.tau_size = len(tau)

        # update min k from class
        return_sources = False

        if sources_array is None:
            # if no array was given a new array is created and returned
            sources_array = np.empty([len(quantities), len(self.get_k()) * self.tau_size], dtype = np.double)
            return_sources = True
        
        # TODO use self.provider?
        start_data_provider = perf_counter()
        provider = CLASSDataProvider(self.cosmo)
        self.times["build DataProvider"] += perf_counter() - start_data_provider #1e-5 sec
        # Get ALL inputs (since we want to be doing this only _once_)

        start_get_inputs = perf_counter()
        raw_inputs = provider.get_inputs(k=self.k, tau=tau, input_selection=self._all_network_input_names()) # 0.002 sec
        transformed_inputs = self.input_transformer.transform_inputs(raw_inputs) # 3e-4 sec

        # Construct a cache object simplifying the access
        cache = PredictorCache(raw_inputs, transformed_inputs) # 3e-6 sec # actually using this cache takes quite a lot of time. Maybe we will find some workaround there ...

        self.times["neural network input transformation"] += perf_counter() - start_get_inputs # 0.014 sec

        start_predict = perf_counter()
        for i in range(len(quantities)):
            sources_array[i] = self.predict(sources_array[i], quantities[i], tau, provider, cache=cache)
        self.times["predictor.predict"] = perf_counter() - start_predict # 0.05 sec


        if return_sources:
            k = self.get_k()
            return k, sources_array
        else:
            return


    def predict_all(self, tau):
        return self.predict_many(["t0", "t1", "t2", "t2_p", "phi_plus_psi", "delta_m", "delta_cb"], tau)

    def untransform(self, quantity, value, raw_inputs):
        start = perf_counter()
        result = self.target_transformer.untransform_target(quantity, value, inputs=raw_inputs)
        elapsed = perf_counter() - start
        self.times["neural network output transformation"] += elapsed

        return result


class ModelWrapper:

    def __init__(self, model):
        self.model = model

    def required_inputs(self):
        raise NotImplementedError

    def __call__(self, inputs):
        pass


class TorchModel(ModelWrapper):

    def __init__(self, model, device, slicing=None):
        self.model = model
        self.device = device
        self.slicing = slicing

    def required_inputs(self):
        return self.model.required_inputs()

    def __call__(self, inputs):
        with torch.no_grad():
            converted_inputs = self._convert_inputs(inputs)
            S_t = self.model(converted_inputs)

        return S_t.cpu().numpy()

    def forward_reduced_mode(self, inputs, k_min_idx):
        with torch.no_grad():
            converted_inputs = self._convert_inputs(inputs)
            S_t = self.model.forward_reduced_mode(converted_inputs, k_min_idx)
        start_trafo = time.time()
        out = S_t.type(torch.DoubleTensor).cpu().numpy()
        trafo_time = time.time()-start_trafo
        return out, trafo_time

    def _convert_inputs(self, inputs):
        # ret = {k: v if v.ndim == 0 else torch.from_numpy(v).float().to(self.device) for k, v in inputs.items() if not isinstance(v, tuple)}
        ret = {k: v if v.ndim == 0 else torch.from_numpy(v.astype(np.float32)).to(self.device) for k, v in inputs.items() if not isinstance(v, tuple)}
        return ret


class TreePredictor(BasePredictor):

    def __init__(self, cosmo, input_transformer, target_transformer, models, rules, funcs=None, k=None):
        super().__init__(cosmo, input_transformer, target_transformer, k=k)

        self.models = models
        self.rules = rules
        self.funcs = funcs if funcs is not None else {}
        self.cache = {}
        self.verbose = False
        self.k_size = len(self.get_k())

    def log(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)

    def update_predictor(self, cosmo, device_name="cpu"):
        """
        If we reuse the predictor instance we need to reload the classy.Class instance to get the updated values.
        """
        # Update cosmo instance
        self.cosmo = cosmo

        # store new minimal value for k, derive new size of output k
        self.k_min_class = self.cosmo.k_min()
        self.k_min_idx = lowest_k_index(self.k, self.k_min_class)
        self.k_size = len(self.get_k())

        # Clean Cache and reset times
        self.reset_cache()
        self.reset_times()

    def _all_network_input_names(self):
        from itertools import chain
        return set(chain(*(mod.required_inputs() for mod in self.models.values())))

    def _predict(self, quantity, tau, provider, cache, source_array):
        cosmo = self.cosmo

        if quantity in self.models:
            source_array = self._predict_from_model(quantity, cosmo, tau, provider, cache, source_array)
        else:
            source_array = self._predict_from_combine(quantity, cosmo, tau, provider, cache, source_array)

        return source_array

    def _predict_from_model(self, quantity, cosmo, tau, provider, cache, source_array, rule = None):
        assert quantity not in self.rules
        if quantity in self.cache:
            start_load_cache = perf_counter()
            S_t, mask_bounds = self.cache[quantity]
            elapsed = perf_counter() - start_load_cache
            self.times["cache load"] += elapsed
        else:
            start_model_load = perf_counter()

            model = self.models[quantity]
            in_select = model.required_inputs()
            elapsed = perf_counter() - start_model_load
            self.times["model load"] += elapsed

            self.log("Evaluating model for", quantity)

            # Check whether we should evaluate model only at certain tau
            slicing = model.slicing
            if slicing is not None:
                mask = slicing.which(cosmo, tau)
                mask_indices = np.where(mask==True)
                mask_bounds = (mask_indices[0][0] * self.k_size,(mask_indices[0][-1] + 1) * self.k_size)
            else:
                mask = None
                mask_bounds = (0,self.tau_size*self.k_size)

            start_input_retrieval = perf_counter()
            inputs = cache.get_transformed_inputs(in_select, tau_mask=mask)

            elapsed = perf_counter() - start_input_retrieval
            self.times["neural network input transformation"] += elapsed

            start_predict = perf_counter()

            #Find lowest index of the NN k-array which is yet to be used. It is until where the input array is to be filled up with
            S_t, trafo_time = model.forward_reduced_mode(inputs, self.k_min_idx) #float type
            elapsed = perf_counter() - start_predict
            self.time_prediction_per_network[quantity] = elapsed - trafo_time
            self.times["float to double transformation"] += trafo_time
            self.times["neural network evaluation"] += elapsed


            #S_t = np.array(S_t,dtype=np.double)
            start_output = perf_counter()

            self.times["neural network output transformation"] += perf_counter() - start_output

            save_cache = perf_counter()
            self.cache[quantity] = (S_t, mask_bounds)
            self.times["save source function in cache"] += perf_counter() - save_cache

        start_output = perf_counter()
        if rule == sum:
            S_t = self._interpolate_k_min(S_t,quantity)
            source_array[mask_bounds[0]:mask_bounds[1]] += np.copy(S_t, order='C')
        elif rule is None:
            S_t = self._interpolate_k_min(S_t,quantity)
            source_array[mask_bounds[0]:mask_bounds[1]] = np.copy(S_t, order='C')
        elif type(rule) == Channel:
            S_t = rule(S_t)
            S_t = self._interpolate_k_min(S_t,quantity)
            source_array[mask_bounds[0]:mask_bounds[1]] = np.copy(S_t, order='C')
        elif rule == 't2_p':
            S_t = self._interpolate_k_min(S_t,quantity)
            source_array[mask_bounds[0]:mask_bounds[1]] += np.copy(S_t*np.sqrt(6), order='C')
        else:
            print("NO RULE FOUND: RIBBIT!")
        self.times["rule application"] += perf_counter() - start_output
        self.times["neural network output transformation"] += perf_counter() - start_output
        
        return source_array

    def _interpolate_k_min(self,S,quantity):
        polyfit_k_min = perf_counter()
        # z = np.polyfit(self.k[self.k_min_idx:self.k_min_idx+3], S[0:3], deg=2)
        # S[0,:] = z[2] + z[1] * self.k_min_class + z[0] * self.k_min_class**2

        # Do log extrapolation here
        # if quantity == ('phi_plus_psi', 'delta_m', 'delta_cb'):
        #     if S[0]<0: #only delta m delta cb is to be interpolated log
        #         m = ( np.log10(-S[2::self.k_size]) - np.log10(-S[1::self.k_size])) / ( np.log10(self.k[self.k_min_idx+2]) - np.log10(self.k[self.k_min_idx+1]) )
        #         S[0::self.k_size] = -10**( np.log10(-S[1::self.k_size]) + m * ( np.log10(self.k_min_class) - np.log10(self.k[self.k_min_idx]) ) )
        #     else:
        #         S[0::self.k_size] = S[1::self.k_size] + (S[2::self.k_size]-S[1::self.k_size])/(self.k[self.k_min_idx+2]-self.k[self.k_min_idx+1]) * (self.k_min_class-self.k[self.k_min_idx])
        # else:
        #     S[0::self.k_size] = S[1::self.k_size] + (S[2::self.k_size]-S[1::self.k_size])/(self.k[self.k_min_idx+2]-self.k[self.k_min_idx+1]) * (self.k_min_class-self.k[self.k_min_idx])
        
        #Do log interpolation here
        # if quantity == ('phi_plus_psi', 'delta_m', 'delta_cb'):

        #     agg = S[:15]
        #     ks = self.k[self.k_min_idx:self.k_min_idx+15]
        #     print(agg)
        #     print(ks)
        #     fig,ax = plt.subplots()
        #     ax.plot(ks,agg)
        #     ax.set_xscale('log')
        #     #ax.set_yscale('log')
        #     fig.savefig('erg.png')
        #     sys.exit()

        #     if S[0]<0: #only delta m delta cb is to be interpolated log
        #         m = ( np.log10(-S[1::self.k_size]) - np.log10(-S[0::self.k_size])) / ( np.log10(self.k[self.k_min_idx+2]) - np.log10(self.k[self.k_min_idx]) )
        #         S[0::self.k_size] = -10**( np.log10(-S[0::self.k_size]) + m * ( np.log10(self.k_min_class) - np.log10(self.k[self.k_min_idx]) ) )
        #     else:
        #         S[0::self.k_size] = S[0::self.k_size] + (S[1::self.k_size]-S[0::self.k_size])/(self.k[self.k_min_idx+1]-self.k[self.k_min_idx+0]) * (self.k_min_class-self.k[self.k_min_idx])
        #     print(S[0::self.k_size])
        #     print(S[1::self.k_size])
        #     print(S[2::self.k_size])
            

        
        
        # else:
        #     S[0::self.k_size] = S[0::self.k_size] + (S[1::self.k_size]-S[0::self.k_size])/(self.k[self.k_min_idx+1]-self.k[self.k_min_idx+0]) * (self.k_min_class-self.k[self.k_min_idx])

        # Do quandratic lengrende interpolation
        k1 = (self.k_min_class-self.k[self.k_min_idx+1]) * (self.k_min_class-self.k[self.k_min_idx+2]) / (self.k[self.k_min_idx] - self.k[self.k_min_idx+1]) / (self.k[self.k_min_idx] - self.k[self.k_min_idx+2])
        k2 = (self.k_min_class-self.k[self.k_min_idx]) * (self.k_min_class-self.k[self.k_min_idx+2]) / (self.k[self.k_min_idx+1] - self.k[self.k_min_idx+2]) / (self.k[self.k_min_idx+1] - self.k[self.k_min_idx])
        k3 = (self.k_min_class-self.k[self.k_min_idx]) * (self.k_min_class-self.k[self.k_min_idx+1]) / (self.k[self.k_min_idx+2] - self.k[self.k_min_idx]) / (self.k[self.k_min_idx+2] - self.k[self.k_min_idx+1])
        S[0::self.k_size] = S[0::self.k_size]*k1 + S[1::self.k_size]*k2 + S[2::self.k_size]*k3


        self.times["interpolate k_min"] += perf_counter() - polyfit_k_min
        return S


    def _predict_from_combine(self, quantity, cosmo, tau, provider, cache, source_array):
        assert quantity not in self.models

        parents, combine, pass_cosmo = self.rules[quantity]
        contributions = []
        self.log("Computing {} as combination of {}.".format(quantity, parents))
        for parent in parents:
            source_array = self._predict_from_model(parent, cosmo, tau, provider, cache, source_array, combine)

        return source_array


def print_dict_sorted(d):
    order = sorted(d, key=lambda key: d[key], reverse=True)
    longest_key = max(len(key) for key in d)
    for key in order:
        print(key.ljust(longest_key + 2), "\t", d[key])

class Timer:
    def __init__(self):
        self._start = {}
        self._times = {}

    def start(self, name):
        self._start[name] = perf_counter()

    def end(self, name):
        self._times[name] = perf_counter() - self._start[name]
        del self._start[name]

    def pprint(self):
        print_dict_sorted(self._times)

# from profilehooks import profile
# @profile(filename="build_predictor.prof")
def build_predictor(cosmo, device_name="cpu"):
    timer = Timer()
    timer.start("all")
    timer.start("cosmo.nn_workspace()")
    workspace = cosmo.nn_workspace()
    timer.end("cosmo.nn_workspace()")
    timer.start("torch.device")
    device = torch.device(device_name)
    timer.end("torch.device")

    timer.start("load k")
    k  = workspace.loader().k()
    timer.end("load k")
    timer.start("move k to device")
    kt = torch.from_numpy(k).float().to(device)
    timer.end("move k to device")
    timer.start("load models")
    models, rules = load_models(workspace, ALL_NETWORK_CLASSES, kt, device)
    timer.end("load models")
    timer.start("build transformer")
    input_transformer, target_transformer = current_transformer.get_pair(workspace.normalization_file, k)
    timer.end("build transformer")

    timer.start("build predictor")
    predictor = TreePredictor(
        cosmo,
        input_transformer, target_transformer,
        models, rules,
        k=k,
    )
    timer.end("build predictor")

    timer.end("all")

    return predictor

def load_models(workspace, classes, k, device):
    from collections import defaultdict
    times = defaultdict(float)
    start_models = perf_counter()
    # models = [ctor(k) for ctor in classes]
    models = []
    for ctor in classes:
        start = perf_counter()
        mod = ctor(k)
        # times[ctor.__name__] = perf_counter() - start
        models.append(mod)
    times["model ctors"] += perf_counter() - start_models
    for model in models:
        start = perf_counter()
        state_dict = torch.load(workspace.model_path(model.name()), map_location=device)
        times["torch.load"] += perf_counter() - start

        start = perf_counter()
        model.load_state_dict(state_dict)
        times["model.load_state_dict"] += perf_counter() - start

        start = perf_counter()
        model.to(device)
        times["model.to(device)"] += perf_counter() - start
        model.eval()

    def model_key(model):
        targets = model.source_functions()
        if len(targets) == 1:
            return targets[0]
        else:
            return tuple(targets)

    def model_wrapper(model):
        return TorchModel(model, device, slicing=model.slicing())

    start = perf_counter()
    model_dict = {model_key(m): model_wrapper(m) for m in models}
    times["model_dict creation"] += perf_counter() - start

    rules = {
            "t0":           (("t0_reco_no_isw", "t0_reio_no_isw", "t0_isw"), sum, False),
            "t2":           (("t2_reco", "t2_reio"), sum, False),
            "t2_p":         (("t2_reco", "t2_reio"), 't2_p', False),
            "phi_plus_psi": ((("phi_plus_psi", "delta_m", "delta_cb"),), Channel(0), False),
            "delta_m":      ((("phi_plus_psi", "delta_m", "delta_cb"),), Channel(1), False),
            "delta_cb":     ((("phi_plus_psi", "delta_m", "delta_cb"),), Channel(2), False),
            }

    return model_dict, rules


class PredictionDescriptor:
    def __init__(self, model_dict, rules):
        """
        Represents how the models need to be evaluated.

        `model_dict` will be a dict from source functions names to `TorchModel`
        instances.

        `rules` will be a dict describing how source functions that are not
        implemented as single networks will be 'assembled' from the outputs
        of other networks.
        For that, `rules` maps the names of the source functions to 3-tuples
        of the form `(tuple_of_dependencies, function_to_apply_to_dependencies, unused)`.
        """
        self.model_dict = model_dict
        self.rules = rules


class Channel:
    def __init__(self, n):
        self.n = n

    def __call__(self, items):
        return items[..., self.n]



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
