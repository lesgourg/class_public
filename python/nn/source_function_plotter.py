import json
from pathlib import Path
import random

import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt

from classy import Class

from .plotting.plot_source_function import plot_source_function

class SourceFunctionPlotter:

    def __init__(self, workspace):
        self.workspace = workspace

    def plot_slice(self, name, marker=None):
        sample, (cosmo, cosmo_nn) = self.get_cosmo_pair()

        sources, k, tau = cosmo.get_sources()
        sources_nn, k_nn, tau_nn = cosmo_nn.get_sources()

        tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
        tau_sel = tau_rec
        i_tau = np.argmin(np.abs(tau - tau_sel))
        i_nn_tau = np.argmin(np.abs(tau_nn - tau_sel))

        fig, ax = plt.subplots()
        ax.semilogx(k, sources[name][:, i_tau], alpha=0.6, label="CLASS", marker=marker)
        ax.semilogx(k_nn, sources_nn[name][:, i_nn_tau], alpha=0.6, label="CLASSnet", marker=marker)
        ax.legend()
        ax.grid()
        ax.set(xlabel="$k$")
        ax.set(ylabel="$S$")
        ax.set_yscale("symlog", linthreshy=0.000001)

        fig.savefig(self.workspace.plots / "slice.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

        cosmo.struct_cleanup()
        cosmo_nn.struct_cleanup()

    def plot_source_functions(self, directory=None):
        # Default to workspace plots directory
        if directory is None:
            directory = self.workspace.plots / "source functions"

        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        assert directory.is_dir()

        sample, (cosmo, cosmo_nn) = self.get_cosmo_pair()

        # Also save the sampled cosmological parameters to the same output directory
        with open(directory / "parameters.json", "w") as out:
            json.dump(sample, out)

        sources, k, tau = cosmo.get_sources()
        sources_nn, k_nn, tau_nn = cosmo_nn.get_sources()
        fields = set(sources) & set(sources_nn)

        tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]

        for field in fields:
            fig = self.create_plot_for_function(
                field,
                k, tau, sources[field],
                k_nn, tau_nn, sources_nn[field],
                tau_rec=tau_rec,
            )
            filename = directory / "{}.png".format(field)
            print("Saving plot of '{}' to '{}'.".format(field, filename))
            fig.savefig(filename, dpi=200, bbox_inches="tight")

        cosmo.struct_cleanup()
        cosmo_nn.struct_cleanup()

    def create_plot_for_function(self, name, k, tau, S, k_nn, tau_nn, S_nn, tau_rec=None):
        fig, axes = plt.subplots(ncols=3, sharey=True, figsize=(3 * 4, 4))
        plot_source_function(
            axes[0],
            k, tau, S,
            tau_rec=tau_rec, levels=50, title="{} (CLASS)".format(name))
        plot_source_function(
            axes[1],
            k_nn, tau_nn, S_nn,
            tau_rec=tau_rec, levels=50, title="{} (CLASSnet)".format(name)
        )

        # compute residual
        spline = RectBivariateSpline(k, tau, S)
        S_on_nn = spline(k_nn, tau_nn)

        residual = S_nn - S_on_nn

        plot_source_function(axes[2], k_nn, tau_nn, residual, tau_rec=tau_rec, levels=50, title="residual")

        return fig

    def get_cosmo_pair(self):
        manifest = self.workspace.loader().manifest()

        _, validation = self.workspace.loader().cosmological_parameters()
        n_validation =  len(validation[next(iter(validation))])
        random_index = random.randrange(0, n_validation)
        sample = {k: v[random_index] for k, v in validation.items()}

        # # TODO REMOVE!!!!
        # val_transpose = [{key: validation[key][i] for key in validation} for i in range(n_validation)]
        # sample = min(val_transpose, key=lambda item: np.abs(item["Omega_k"] - (-0.011)))
        # print("SAMPLE:", sample)

        params = {}
        params.update(manifest["fixed"])
        params.update(sample)

        cosmo = Class()
        cosmo.set(params)
        cosmo.compute(level=["perturb"])

        cosmo_nn = Class()
        cosmo_nn.set(params)
        cosmo_nn.set({"neural network path": self.workspace, "nn_debug": True})
        cosmo_nn.compute(level=["perturb"])

        return sample, (cosmo, cosmo_nn)

