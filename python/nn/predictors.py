import multiprocessing
import functools

from .predictor_cache import PredictorCache

from .data_providers import CLASSDataProvider
# from plotting.plot_source_function import plot_source_function
# import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import os
from time import perf_counter
import scipy.interpolate

import torch

from .models import ALL_NETWORK_CLASSES
from . import current_transformer

def lowest_k_index(k, k_min):
    idx = np.searchsorted(k, k_min)
    return max(idx - 1, 0)

class BasePredictor:
    def __init__(self, cosmo, input_transformer, target_transformer, k):
        self.cosmo = cosmo

        self.k = k

        self.input_transformer = input_transformer
        self.target_transformer = target_transformer

        self.reset_cache()
        self.reset_times()

    def reset_cache(self):
        self.raw_cache = {}
        self.transformed_cache = {}

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
        else:
            return 0

    def predict(self, quantity, tau, provider=None, cache=None):
        """
        Get the network prediction for the given `quantity` and `tau` array.
        The result will be sampled on the `k` array obtained from `self.get_k()` and
        thus have the shape `(len(self.get_k()), len(tau))`
        """
        start_predict = perf_counter()
        result, raw_inputs = self._predict(quantity, tau, provider, cache)
        k_min_class = self.cosmo.k_min()

        # NOTE: THIS IS VERY IMPORTANT! WE MAY NOT PASS K MODES < K_MIN_CLASS TO CLASS
        # OTHERWISE BAD THINGS WILL HAPPEN!!!
        k_min_idx = lowest_k_index(self.k, k_min_class)
        right = result[self.k > k_min_class]

        # t = (k_min_class - self.k[k_min_idx]) / (self.k[k_min_idx + 1] - self.k[k_min_idx])
        # left = result[k_min_idx] * (1.0 - t) + result[k_min_idx + 1] * t

        z = np.polyfit(self.k[k_min_idx:k_min_idx+3], result[k_min_idx:k_min_idx+3], deg=2)
        left = z[2] + z[1] * k_min_class + z[0] * k_min_class**2

        # old_result = result
        result = np.concatenate((left[None, :], right), axis=0)

        k = self.get_k()

        # here we handle the special case of delta_m by fitting the (k s_2)^2 factor
        # for low values of k, i.e. delta_m = p * (k s_2)^2

        if quantity == "delta_m":
            S = result
            assert S.shape == (len(k), len(tau))
            print("INFO: delta_m")
            k_th = 5e-4
            fit_mask = (k > k_th) & (k < 1e-3)
            from scipy.optimize import leastsq
            Omega_k = raw_inputs["cosmos/Omega_k"]
            h = raw_inputs["cosmos/h"]
            def get_factor(k):
                return k[:, None]**2 + 3 * Omega_k * (h / 2997.9)**2
            slope = S[fit_mask].sum(axis=0) / get_factor(k[fit_mask]).sum(axis=0)
            # slope = S[fit_mask][0] / get_factor(k[fit_mask])[0]
            replacement = slope * get_factor(k[k < k_th])

            # import matplotlib; matplotlib.use("qt4agg")
            # import matplotlib.pyplot as plt
            # plt.loglog(k, -S[:, -1], label="prediction", color="g")
            # plt.loglog(k[k < k_th], -replacement[:, -1], label="replacement", ls="-.")
            # plt.axvspan(k[fit_mask][0], k[fit_mask][-1], color="yellow", alpha=0.4)
            # # plt.loglog(k < k_th], -replacement[:, -1], label="replacement", ls="-.", c="r")
            # plt.legend()
            # plt.grid()
            # plt.show()

            result[k < k_th] = replacement


        # i_tau = 120
        # import matplotlib
        # matplotlib.use("qt5agg")
        # import matplotlib.pyplot as plt
        # plt.semilogx(self.k, old_result[:, i_tau], label="network prediction")
        # plt.semilogx(self.k[self.k > k_min_class], right[:, i_tau], label="right part", ls="--")
        # plt.scatter([k_min_class], [left[i_tau]], label="point at k_min_class")
        # plt.semilogx(self.get_k(), result[:, i_tau], label="final result", ls="-.")
        # plt.axvline(k_min_class, label="k_min_class", color="r", ls="--")
        # plt.axvline(self.k[k_min_idx], label="idx one before k_min_class", color="g", ls="--")
        # plt.scatter(self.k[k_min_idx:k_min_idx+3], old_result[k_min_idx:k_min_idx+3, i_tau], label="points for quadratic")
        # plt.legend()
        # plt.grid()
        # plt.show()

        self.times["predictor.predict"] += perf_counter() - start_predict
        return result
        # return result[lowest_k_index(self.k, k_min_class), :]

    def _predict(self, quantity, tau):
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
        k_min_class = self.cosmo.k_min()

        # NOTE: THIS IS VERY IMPORTANT! WE MAY NOT PASS K MODES < K_MIN_CLASS TO CLASS
        # OTHERWISE BAD THINGS WILL HAPPEN!!!
        k_result = np.concatenate((
            np.array([k_min_class]),
            k[k > k_min_class]
        ))
        return k_result

    def predict_many(self, quantities, tau):
        """
        Predict the source functions whose names are given as the list `quantities`.
        This will return a numpy array of shape (len(quantities), len(k) + 1, len(tau)).
        The 2nd size is len(k) + 1 instead of len(k) because this function adds another
        row to the S array corresponding to k -> 0.
        This is needed so that CLASS does not have to extrapolate for low k.
        """
        # TODO use self.provider?
        provider = CLASSDataProvider(self.cosmo)

        start = perf_counter()
        # Get ALL inputs (since we want to be doing this only _once_)
        raw_inputs = provider.get_inputs(k=self.k, tau=tau, input_selection=self._all_network_input_names())
        transformed_inputs = self.input_transformer.transform_inputs(raw_inputs)

        # Construct a cache object simplifying the access
        cache = PredictorCache(raw_inputs, transformed_inputs)

        self.times["neural network input transformation"] += perf_counter() - start

        predictions = {qty: self.predict(qty, tau, provider, cache=cache) for qty in quantities}

        k = self.get_k()
        k_len = len(k)
        result = np.zeros((len(quantities), len(k), len(tau)))

        # Store predictions in array
        for i, quantity in enumerate(quantities):
            S = predictions[quantity]
            result[i, :, :] = S

        return k, result

    def predict_all(self, tau):
        return self.predict_many(["t0", "t1", "t2", "phi_plus_psi", "delta_m"], tau)

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
        # self.model.eval()
        from time import perf_counter
        with torch.no_grad():
            a = perf_counter()
            converted_inputs = self._convert_inputs(inputs)
            b = perf_counter()
            S_t = self.model(converted_inputs)
            c = perf_counter()

            t_convert = b - a
            t_infer = c - b
            t_tot = c - a

            # print("convert:\t{}s\ninfer:\t{}s\ntotal:\t{}s".format(
            #     t_convert,
            #     t_infer,
            #     t_tot
            #     ))

        return S_t.cpu().numpy()

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

    def log(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)

    def _all_network_input_names(self):
        from itertools import chain
        return set(chain(*(mod.required_inputs() for mod in self.models.values())))

    def _predict(self, quantity, tau, provider, cache):
        cosmo = self.cosmo

        if quantity in self.models:
            S, raw_inputs = self._predict_from_model(quantity, cosmo, tau, provider, cache)
        else:
            S, raw_inputs = self._predict_from_combine(quantity, cosmo, tau, provider, cache)

        # TODO CAREFUL
        if quantity == "t2_reco":
            print("T2RECO EXTRAPOLATION")
            k_th = 1e-3
            i_k_th = np.argmin(np.abs(self.k - k_th))
            S_th = S[i_k_th]
            k_mask = self.k < k_th
            S[k_mask] = S_th * (self.k[k_mask] / k_th)[:, None]**2


        if cosmo.nn_cheat_enabled() and quantity in cosmo.nn_cheat_sources():
            # Here, the `quantity` should not be predicted by a network, but
            # instead be taken from CLASS.
            # First, emit a warning message to make sure that this is not
            # accidentally enabled
            print("WARNING: 'CHEATING' IS ENABLED FOR QUANTITY '{}'".format(quantity))
            # It is guaranteed that if `cosmo.nn_cheat_mode()` is true,
            # the perturbation module has been fully executed and we can
            # thus simply take the source function in question from there.
            # This is done after evaluating a networks, because we also need
            # to return the raw inputs, i.e. we only replace the source function.
            S_cheat, k_cheat, tau_cheat = cosmo.get_sources()
            # it remains to perform interpolation of the source function onto
            # the desired (k, tau)-grid.

            spline = scipy.interpolate.RectBivariateSpline(k_cheat, tau_cheat, S_cheat[quantity])
            # this function must also return the `raw_inputs` dict
            return spline(self.k, tau), raw_inputs

            # TODO TODO TODO
            # k_th = 1e-10
            # k_left = k_cheat[k_cheat < k_th]
            # k_right = self.k[self.k >= k_th]
            # k_ = np.concatenate((k_left, k_right))
            # S_left = S_cheat[quantity][k_cheat < k_th]
            # S_right = S[self.k >= k_th]
            # S_ = np.concatenate((S_left, S_right))
            # spline = scipy.interpolate.RectBivariateSpline(k_, tau, S_)
            # # this function must also return the `raw_inputs` dict
            # return spline(self.k, tau), raw_inputs

        if quantity in self.funcs:
            self.funcs[quantity](S, raw_inputs)

        return S, raw_inputs

    def _predict_from_model(self, quantity, cosmo, tau, provider, cache):
        assert quantity not in self.rules
        if quantity in self.cache:
            return self.cache[quantity]

        model = self.models[quantity]

        in_select = model.required_inputs()

        self.log("Evaluating model for", quantity)

        # Check whether we should evaluate model only at certain tau
        slicing = model.slicing
        if slicing is not None:
            mask = slicing.which(cosmo, tau)
            tau_eval = tau[mask]
        else:
            tau_eval = tau
            mask = None

        start_input_retrieval = perf_counter()
        raw_inputs = cache.get_raw_inputs(in_select, tau_mask=mask)
        inputs = cache.get_transformed_inputs(in_select, tau_mask=mask)
        elapsed = perf_counter() - start_input_retrieval
        self.times["neural network input transformation"] += elapsed

        start_predict = perf_counter()
        S_t = model(inputs)
        elapsed = perf_counter() - start_predict
        self.time_prediction_per_network[quantity] = elapsed
        self.times["neural network evaluation"] += elapsed

        result_shape = list(S_t.shape)
        result_shape[0] = tau.size
        result = np.zeros(result_shape)

        if slicing is not None:
            result[mask] = S_t
            self.log("Slicing output for quantity", quantity)
        else:
            result[:, :] = S_t

        S_t = result


        start_output = perf_counter()
        # Swap k and tau axis
        S_t = np.swapaxes(S_t, 0, 1)
        if isinstance(quantity, tuple) or isinstance(quantity, list):
            # If quantity is a 'container' for multiple quantities (e.g. phi+psi and delta_m),
            # transform the components individually
            S = np.stack([self.untransform(q, S_t[..., i], raw_inputs) for i, q in enumerate(quantity)], axis=2)
        else:
            S = self.untransform(quantity, S_t, raw_inputs)
        self.times["neural network output transformation"] += perf_counter() - start_output

        self.cache[quantity] = (S, raw_inputs)
        return S, raw_inputs

    def _predict_from_combine(self, quantity, cosmo, tau, provider, cache):
        assert quantity not in self.models
        parents, combine, pass_cosmo = self.rules[quantity]

        contributions = []
        raw_inputs = {}
        self.log("Computing {} as combination of {}.".format(quantity, parents))
        for parent in parents:
            contrib, raw_inp = self._predict(parent, tau, provider, cache)
            contributions.append(contrib)
            raw_inputs.update(raw_inp)

        if pass_cosmo:
            S = combine(contributions, cosmo, tau)
        else:
            S = combine(contributions)

        return S, raw_inputs


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
    timer.pprint()

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

    print("load_models:")
    print_dict_sorted(times)

    rules = {
            "t0":           (("t0_reco_no_isw", "t0_reio_no_isw", "t0_isw"), sum, False),
            "t2":           (("t2_reco", "t2_reio"), sum, False),
            "phi_plus_psi": ((("phi_plus_psi", "delta_m"),), Channel(0), False),
            "delta_m":      ((("phi_plus_psi", "delta_m"),), Channel(1), False),
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
        return items[0][..., self.n]
