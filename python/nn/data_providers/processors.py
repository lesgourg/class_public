import numpy as np
import scipy.interpolate

from .containers import InputContainer, TargetContainer
from .approximations import add_approximations

import utils

class InputPostProcessor:

    def process_default(self, container, k, tau=None):
        self.process_scalars(container)

        # TODO Is this necessary or is it more useful to treat this at the level of the network?
        ones = np.ones_like(tau)
        scalars = {k: ones * v for k, v in container.scalars.items()}
        cosmos = {k: ones * v for k, v in container.cosmos.items()}

        background = self.process_background(container, tau)
        thermo = self.process_thermo(container, tau)

        result = dict()
        result.update(cosmos)
        result.update(scalars)
        result.update(background)
        result.update(thermo)

        result["k"] = np.tile(k[None, :], (len(tau), 1))
        result["tau"] = tau

        self.add_relative_times(container, tau, result)

        return result

    def process(self, container, k, tau=None, selection=None):
        assert isinstance(container, InputContainer)
        tau = tau if tau is not None else container.tau_source

        result = self.process_default(container, k, tau)

        add_approximations(result, k, tau, selection=selection)

        return result

    def add_relative_times(self, container, tau, result):
        result["tau_relative_to_reco"] = tau / container.scalars.tau_rec
        result["tau_relative_to_reio"] = tau / container.scalars.tau_reio

    def process_scalars(self, container):
        # container.scalars.k_eq = container.scalars.a_eq * container.scalars.H_eq

        container.scalars.tau_eisw_lisw_50_relative_to_rec = container.scalars.tau_eisw_lisw_50 / container.scalars.tau_rec
        container.scalars.tau_eisw_lisw_120_relative_to_rec = container.scalars.tau_eisw_lisw_120 / container.scalars.tau_rec

        container.scalars.tau_eisw_lisw_50_relative_to_reio = container.scalars.tau_eisw_lisw_50 / container.scalars.tau_reio
        container.scalars.tau_eisw_lisw_120_relative_to_reio = container.scalars.tau_eisw_lisw_120 / container.scalars.tau_reio

    def process_background(self, container, tau):
        background = dict()
        for key, value in container.background.items():
            background[key] = np.interp(tau, container.tau_bg, value)
        return background

    def process_thermo(self, container, tau):
        thermo = dict()
        tau_th_flipped = np.flip(container.tau_th)
        for key, value in container.thermo.items():
            thermo[key] = np.interp(tau, tau_th_flipped, np.flip(value))
        if "r_d" in thermo:
            thermo["k_d"] = 2 * np.pi / thermo["r_d"]
            del thermo["r_d"]
        return thermo

class TargetPostProcessor:
    def __init__(self):
        pass

    def process(self, container, k, tau=None):
        assert isinstance(container, TargetContainer)

        tau = tau if tau is not None else container.tau_source

        result = dict()

        for key, S in container.sources.items():
            # perform extrapolation (very important for delta_m)
            # this is necessary since some source functions are sampled
            # on k arrays that are true subsets of `k_standard`; feeding
            # those functions into RectBivariateSpline will lead to nearest
            # neighbor interpolation, causing our networks to learn wrong
            # values outside the source functions' sampling domain.
            # k_source, S = utils.extrapolate_pow(container.k_source, S.T, k)
            # EDIT: No longer necessary; training data is now always generated
            # such that extrapolation never has to occur
            spline = scipy.interpolate.RectBivariateSpline(container.tau_source, container.k_source, S.T)
            S_interpolated = spline(tau, k)
            result[key] = S_interpolated

        return result
