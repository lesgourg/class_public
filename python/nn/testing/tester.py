import multiprocessing
import random

import numpy as np
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt

from tqdm import tqdm
import classy
from classy import Class

from ..generate.generate_cosmological_parameters import sample_cosmological_parameters

# style for error plots
LINESTYLE = dict(ls="-", lw=0.5, alpha=0.4, color="r")

class Tester:
    def __init__(self, workspace):
        self.workspace = workspace

        self.k = self.workspace.loader().k()
        manifest = workspace.loader().manifest()

        self.fixed = manifest["fixed"]

        self.init_plots()

    def test(self, count=None, cheat=None, prefix=None):
        """
        test the performance of the trained networks by computing the observables
        (Cl's and P(k)) for `count` cosmological scenarios.
        For each observable a plot will be created that displays the deviation
        of the network prediction from the true value of the observable.

        Optionally, `cheat` can be provided as a list of names of source functions
        for which to ignore the network predictions and instead take the true
        value from CLASS.
        This aids in isolating networks which may not perform well.
        """

        _, params = self.workspace.loader().cosmological_parameters()
        validation_size = len(params[next(iter(params))])
        if count is None:
            count = validation_size

        if count > validation_size:
            print("WARNING: You requested a validation same count of {0}, \
                  but there are only {1}! Proceeding with {1}".format(count, validation_size))
            count = validation_size

        selection = random.sample(range(validation_size), count)
        print("Evaluating network predictions for {} cosmologies.".format(count))

        cosmo_params = [{k: v[i] for k, v in params.items()} for i in selection]
        # params = sample_cosmological_parameters(self.domain, count)

        class_params = list(map(self.get_class_parameters, cosmo_params))
        class_params_nn = list(map(self.get_class_parameters_nn, cosmo_params))

        # for recording cl/pk errors as functions of cosmological parameters
        stats = []

        class_pair_iter = map(self.run_class_for, zip(class_params, class_params_nn))
        counter = 0
        for cosmo_param_dict, (exc, pair) in zip(cosmo_params, class_pair_iter):
            if exc:
                print("An exception occured:")
                print(str(exc))
                print("continuing...")
                continue
            (cl_true, k_pk, pk_true), (cl_nn, _, pk_nn) = pair
            counter += 1
            print("PROGRESS: {}".format(counter))

            self.update_plots(cl_true, cl_nn, k_pk, pk_true, pk_nn)
            self.update_stats(stats, cosmo_param_dict, cl_true, cl_nn, k_pk, pk_true, pk_nn)

        self.save_stats(stats, prefix=prefix)
        self.save_plots(prefix=prefix)

    def update_stats(self, stats, cosmo_params, cl_true, cl_nn, k_pk, pk_true, pk_nn):
        cl_err = {key: cl_nn[key] - cl_true[key] for key in cl_true}
        pk_err = pk_nn - pk_true

        cl_err_relative = {key: (cl_nn[key] - cl_true[key]) / cl_true[key] for key in cl_true}
        pk_err_relative = (pk_nn - pk_true) / pk_true

        stat_dict = {
            "parameters": cosmo_params,
            "cl_error": cl_err,
            "k_pk": k_pk,
            "pk_error": pk_err,
            "cl_error_relative": cl_err_relative,
            "pk_error_relative": pk_err_relative
        }

        # stat_dict = {k: v if not isinstance(v, np.ndarray) else list(v) for k, v in stat_dict.items()}
        stats.append(stat_dict)

    def get_class_parameters(self, cosmo_params):
        """ add fixed parameters (no nn) """
        params = {}
        params.update(cosmo_params)
        params.update(self.fixed)
        return params

    def get_class_parameters_nn(self, cosmo_params, cheat=None):
        """ add fixed parameters (with nn) """
        params = self.get_class_parameters(cosmo_params)
        params["neural network path"] = self.workspace
        # if cheating is enabled, we need to notify classy of the quantities
        if cheat:
            params["nn_cheat"] = cheat
        return params

    def parameter_dicts(self, cosmo_params, cheat=None):
        """
        iterator over tuples of `(params_without_nn, params_with_nn)` for
        given domain.
        """
        count = len(cosmo_params)
        for i in range(count):
            params = {}
            params.update(self.fixed)
            params.update({k: v[i] for k, v in cosmo_params.items()})

            params_nn = params.copy()
            params_nn["neural network path"] = self.workspace
            # if cheating is enabled, we need to notify classy of the quantities
            if cheat:
                params_nn["nn_cheat"] = cheat

            yield params, params_nn

    def run_class_for(self, param_tuple):
        """
        run class twice: once without nn, and once with.
        """
        params, params_nn = param_tuple
        try:
            ret = self.run_class(params), self.run_class(params_nn)
            exc = None
        except classy.CosmoComputationError as e:
            ret = None
            exc = e

        return exc, ret

    def k_pk(self, cosmo):
        k_min = cosmo.k_min()
        k = self.k[self.k > k_min]
        k = np.insert(k, 0, k_min)
        return k

    def run_class(self, params):
        """
        Run CLASS for the given `params` and return the Cls and mPk.
        """
        cosmo = Class()
        cosmo.set(params)
        cosmo.compute()
        cls = truncate(cosmo.raw_cl())
        k_pk = self.k_pk(cosmo)
        # TODO maybe insert this check also into cosmo.pk?
        assert np.all(k_pk >= cosmo.k_min())
        pk = np.vectorize(cosmo.pk)(k_pk, 0.0)
        cosmo.struct_cleanup()

        return cls, k_pk, pk

    def init_plots(self):
        self.figs = {
            "tt": plt.subplots(),
            "ee": plt.subplots(),
            "te": plt.subplots(),
            "pk": plt.subplots()
        }

        for _, ax in self.figs.values():
            ax.grid()

        for q in ("tt", "te", "ee"):
            fig, ax = self.figs[q]
            fig.tight_layout()
            ax.set_xlabel(r"$\ell$")

        self.figs["tt"][1].set_ylabel(r"$\Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        ll1 = r"$\ell (\ell + 1) "
        self.figs["te"][1].set_ylabel(ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["ee"][1].set_ylabel(ll1 + r"\Delta C_\ell^{EE}$")

        self.figs["pk"][0].tight_layout()
        self.figs["pk"][1].set(
            xlabel=r"$k$ [Mpc${}^{-1}$]",
            ylabel=r"$(P_{NN}(k)-P_{CLASS}(k))/P_{CLASS}(k)$",
        )

    def update_plots(self, cl_true, cl_nn, k_pk, pk_true, pk_nn):
        ell = cl_true["ell"]

        for _, ax in self.figs.values():
            ax.axhline(0, color="k")

        # TT
        tt_relerr = (cl_nn["tt"] - cl_true["tt"]) / cl_true["tt"]
        self.figs["tt"][1].semilogx(ell, tt_relerr, **LINESTYLE)

        # TE + EE
        for q in ("ee", "te"):
            err = (cl_nn[q] - cl_true[q])
            err_relmax = err / cl_true[q].max()
            self.figs[q][1].semilogx(ell, ell * (ell + 1) * err, **LINESTYLE)

        # P(k)
        pk_relerr = (pk_nn - pk_true) / pk_true
        self.figs["pk"][1].semilogx(k_pk, pk_relerr, **LINESTYLE)

    def save_stats(self, stats, prefix=None):
        import pickle
        fname =  "errors.pickle"
        if not prefix:
            path = self.workspace.plots / fname
        else:
            dir_path = self.workspace.plots / prefix
            dir_path.mkdir(parents=True, exist_ok=True)
            path = dir_path / fname
        print("Writing Cl/Pk error stats to '{}'".format(path))
        with open(path, "wb") as out:
            pickle.dump(stats, out)

    def save_plots(self, prefix=None):
        for name, (fig, _) in self.figs.items():
            if not prefix:
                path = self.workspace.plots / "{}.png".format(name)
            else:
                dir_path = self.workspace.plots / prefix
                dir_path.mkdir(parents=True, exist_ok=True)
                path = dir_path / "{}.png".format(name)

            print("saving plot to", path)
            fig.savefig(path, dpi=200, bbox_inches="tight")

def truncate(cls):
    """
    Removes the first two multipoles of a Cl dict.
    """
    return {k: v[2:] for k, v in cls.items()}
