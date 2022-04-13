import numpy as np
import scipy.interpolate
from .containers import InputContainer, TargetContainer
from .approximations import add_approximations

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
            # Do a quadratic extrapolation using the lowest 3 points of the source function
            # in k.
            deg = 2
            z = np.polyfit(container.k_source[:deg+1], S[:deg+1], deg=deg)
            k_np = k.numpy()
            # a bit of tolerance because otherwise two k points will
            # be extremely close leading to numerical problems in the
            # spline further below
            rtol = 1e-3
            k_low = k_np[k_np < (1.0 - rtol) * container.k_source[0]]
            S_low = z[2] + z[1] * k_low[:, None] + z[0] * k_low[:, None]**2
            k_interp = np.concatenate((k_low, container.k_source), axis=0)
            S_new    = np.concatenate((S_low, S), axis=0)
            assert S_new.shape == (len(k_interp), len(tau))

            # import matplotlib as mpl
            # mpl.use("qt5agg")
            # import matplotlib.pyplot as plt
            # plt.figure()
            # print("source function:", key)
            # index = 118
            # plt.semilogx(container.k_source, S[:, index], label="from CLASS")
            # plt.semilogx(k_low, S_low[:, index], marker="x", label="extrapolated")
            # plt.scatter(container.k_source[:deg+1], S[:deg+1, index], color="r", label="points used in extraploation")
            # plt.legend()
            # plt.grid()
            # plt.show()

            # finally, perform interpolation onto `k`
            spline = scipy.interpolate.RectBivariateSpline(container.tau_source, k_interp, S_new.T)
            S_interpolated = spline(tau, k)
            result[key] = S_interpolated

        return result
