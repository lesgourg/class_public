import multiprocessing

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

    def test(self, count, processes=None, cheat=None, prefix=None):
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

        print("Evaluating network predictions for {} cosmologies.".format(count))
        _, params = self.workspace.loader().cosmological_parameters()
        params = {k: v[:count] for k, v in params.items()}
        # params = sample_cosmological_parameters(self.domain, count)

        # with multiprocessing.Pool(processes) as pool:
        #     iter_ = pool.imap_unordered(self.run_class_for, self.parameter_dicts(params, cheat=cheat))
        iter_ = map(self.run_class_for, self.parameter_dicts(params, cheat=cheat))
        for _ in range(1):
            for exc, pair in tqdm(iter_, total=count):
                if exc:
                    print("An exception occured:")
                    print(str(exc))
                    print("continuing...")
                    continue
                (cl_true, pk_true), (cl_nn, pk_nn) = pair

                self.update_plots(cl_true, cl_nn, pk_true, pk_nn)

        self.save_plots(prefix=prefix)

    def parameter_dicts(self, cosmo_params, cheat=None):
        """
        iterator over tuples of `(params_without_nn, params_with_nn)` for
        given domain.
        """
        count = len(cosmo_params[next(iter(cosmo_params))])
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

    def run_class(self, params):
        """
        Run CLASS for the given `params` and return the Cls and mPk.
        """
        cosmo = Class()
        cosmo.set(params)
        cosmo.compute()
        cls = truncate(cosmo.raw_cl())
        pk = np.vectorize(cosmo.pk)(self.k, 0.0)
        cosmo.struct_cleanup()

        return cls, pk

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

    def update_plots(self, cl_true, cl_nn, pk_true, pk_nn):
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
        self.figs["pk"][1].semilogx(self.k, pk_relerr, **LINESTYLE)

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
