import argparse
import pickle
from pathlib import Path
from abc import ABC, abstractmethod

import numpy as np

# TODO should probably not do this globally
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

PRETTY_LABELS = {
    "N_ur":       r"$N_{ur}$",
    "Omega_k":    r"$\Omega_k$",
    "h":          r"$h$",
    "omega_b":    r"$\omega_b$",
    "omega_cdm":  r"$\omega_\mathrm{cdm}$",
    "omega_ncdm": r"$\omega_\mathrm{ncdm}$",
    "tau_reio":   r"$\tau_\mathrm{reio}$",
    "w0_fld":     r"$w_0$",
    "wa_fld":     r"$w_a$",
}

def get_pretty_label(name):
    return PRETTY_LABELS.get(name, name)

def rms(x):
    return np.sqrt(np.mean(np.square(x)))

class ErrorQuantifier(ABC):

    def __init__(self, title):
        self.title = title

    @abstractmethod
    def __call__(self, item):
        pass

    @abstractmethod
    def file_name(self):
        pass

class RelativeRMSErrorCl(ErrorQuantifier):

    def __init__(self, title, field):
        super().__init__(title)
        self.field = field

    def __call__(self, item):
        return rms(item["cl_error_relative"][self.field])

    def file_name(self):
        return "cl_{}_rel_err_rms".format(self.field)


class RelativeRMSErrorPk(ErrorQuantifier):

    def __call__(self, item):
        return rms(item["pk_error_relative"])

    def file_name(self):
        return "pk_rel_err_rms"


class RelativeMaxErrorCl(ErrorQuantifier):

    def __init__(self, title, field):
        super().__init__(title)
        self.field = field

    def __call__(self, item):
        return np.max(np.abs(item["cl_error_relative"][self.field]))

    def file_name(self):
        return "cl_{}_rel_err_max".format(self.field)


class RelativeMaxErrorPk(ErrorQuantifier):

    def __call__(self, item):
        return np.max(np.abs(item["pk_error_relative"]))

    def file_name(self):
        return "pk_rel_err_max"

class TrianglePlotter:
    def __init__(self, workspace):
        self.workspace = workspace

        self.quantifiers = [
            RelativeRMSErrorCl("RMS relative error of $C_\ell^{TT}$ (log10)", "tt"),
            # TODO ^ relative to CV
            RelativeMaxErrorCl("max. relative error (over all $\\ell$) of $C_\\ell^{TT}$ (log10)", "tt"),
            # TODO ^ relative to CV
            RelativeRMSErrorPk("RMS of relative error (over all $k$) of $P(k)$ (log10)"),
            RelativeMaxErrorPk("max. relative error (over all $k$) of $P(k)$ (log10)")
        ]

    def plot_and_save(self):
        stats = self.workspace.loader().stats()

        param_names = list(stats[0]["parameters"].keys())
        param_values = np.array(
            [[item["parameters"][name] for name in param_names] for item in stats]
        )

        print("creating scatter plots of "
              "cosmological parameters vs. spectra errors")

        for quantifier in self.quantifiers:
            errors = np.array([quantifier(item) for item in stats])
            errors = np.log10(errors)
            fig = self.plot(param_names, param_values, errors)
            fig.suptitle(quantifier.title)
            save_path = self.workspace.plots / (quantifier.file_name() + ".pdf")
            print("saving triangle plot for '{}' to {}".format(quantifier.title, save_path))
            fig.savefig(save_path)
            plt.close(fig)

    def plot(self, param_names, param_values, errors, rasterized=True):
        """
        param_names: List[str]
        param_values: List[List[float]] <- inner list must have same order as param_names
        errors: List[float]
        """

        size = len(param_names)
        fig, axes = plt.subplots(nrows=size, ncols=size, figsize=(3 * size, 3 * size))
        fig.subplots_adjust(wspace=0, hspace=0)

        err_mean = np.mean(errors)
        err_std = np.std(errors)

        err_color = errors

        # remove plots above diagonal
        for i, j in zip(*np.triu_indices_from(axes, 1)):
            axes[i, j].set_visible(False)

        for i, j in zip(*np.tril_indices_from(axes)):
            ax = axes[i, j]
            if i == size - 1:
                param_x = get_pretty_label(param_names[j])
                ax.set_xlabel(param_x)
            if j == 0:
                param_y = get_pretty_label(param_names[i])
                ax.set_ylabel(param_y)
            ax.label_outer()

            if i == j:
                x = param_values[:, j]
                ax_t = ax.twinx()
                ax_t.scatter(x, errors, c=err_color, rasterized=rasterized)
                ax.set_xlim(x.min(), x.max())
                if i == 0:
                    ax.set_ylim(x.min(), x.max())
                else:
                    ax.get_yaxis().set_ticks([])
            else:
                x = param_values[:, j]
                y = param_values[:, i]
                ax.scatter(x, y, c=err_color, rasterized=rasterized)
                ax.set_xlim(x.min(), x.max())
                ax.set_ylim(y.min(), y.max())

        return fig
