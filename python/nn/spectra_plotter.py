import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
import matplotlib.pyplot as plt

# style for error plots
LINESTYLE_RED = dict(ls="-", lw=0.5, alpha=0.4, color="r")
LINESTYLE_GREEN = dict(ls="-", lw=0.5, alpha=0.4, color="g")
LINESTYLE_BLUE = dict(ls="-", lw=0.5, alpha=0.4, color="b")

class SpectraPlotter:

    def __init__(self, workspace):
        self.workspace = workspace

    def plot(self, include_params=False):
        print("plotting spectra...")
        self._init()
        stats = self.workspace.loader().stats()
        print("building plots...")

        # this is terribly slow but whatever...
        for row in stats:
            self._update_plots(
                row["cl_true"], row["cl_nn"],
                row["k_pk"], row["pk_true"],
                row["k_pk_nn"], row["pk_nn"]
            )

        if include_params:
            print("creating color plot for TT")
            self._plot_colored(stats, field="tt")
            print("creating color plot for P(k)")
            self._plot_colored(stats, field="pk")

        self._save()
        self._close_figs()

    def _init(self):
        self.figs = {
            "tt": plt.subplots(),
            # cosmic variance
            "tt_cv": plt.subplots(),
            "ee": plt.subplots(),
            "te": plt.subplots(),
            "pk": plt.subplots(),
            "pk_abs": plt.subplots()
        }

        for _, ax in self.figs.values():
            ax.grid()

        for q in ("tt", "tt_cv", "te", "ee"):
            fig, ax = self.figs[q]
            fig.tight_layout()
            ax.set_xlabel(r"$\ell$")

        self.figs["tt"][1].set_ylabel(r"$\Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        self.figs["tt_cv"][1].set_ylabel(r"$\sqrt{\frac{2\mathrm{min}(\ell, 2000)+1}{2}} \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        ll1 = r"$\ell (\ell + 1) "
        self.figs["te"][1].set_ylabel(ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["ee"][1].set_ylabel(ll1 + r"\Delta C_\ell^{EE}$")
        self.figs["ee"][1].set_yscale("symlog", linthresh=1e-14)

        self.figs["pk"][0].tight_layout()
        self.figs["pk"][1].set(
            xlabel=r"$k$ [Mpc${}^{-1}$]",
            ylabel=r"$(P_{NN}(k)-P_{CLASS}(k))/P_{CLASS}(k)$",
            xlim=(1e-4, 100),
        )

        self.figs["pk_abs"][0].tight_layout()
        self.figs["pk_abs"][1].set(
            xlabel=r"$k$ [Mpc${}^{-1}$]",
            ylabel=r"$P_{NN}(k)$",
        )

    def _get_param_bounds(self, data):
        param_names = data[0]["parameters"].keys()
        params = {key: np.array([item["parameters"][key] for item in data]) for key in param_names}
        return {key: (np.min(val), np.max(val)) for key, val in params.items()}

    def _plot_colored(self, data, field="tt"):
        def plot_item_tt(ax, item, **kwargs):
            ax.set_xlabel("$\\ell$")
            ax.set_ylabel("$\\Delta C_\\ell^\mathrm{TT} / C_\\ell^\mathrm{TT}$")
            ax.semilogx(item["cl_true"]["ell"], item["cl_error_relative"]["tt"], **kwargs)

        def plot_item_pk(ax, item, **kwargs):
            ax.set_xlabel("$k$ [Mpc${}^{-1}]$")
            ax.set_ylabel("$\\Delta P(k) / P(k)$")
            ax.semilogx(item["k_pk"], item["pk_error_relative"], **kwargs)
            ax.set_yscale("symlog", linthresh=0.01)

        handlers = {
            "tt": plot_item_tt,
            "pk": plot_item_pk,
        }
        if field not in handlers:
            raise ValueError("invalid field: `{}`".format(tt))
        else:
            handler = handlers[field]

        pbounds = self._get_param_bounds(data)
        cmap = plt.cm.winter
        for param_name, bound in pbounds.items():
            print("_plot_cl_colored({}) for {}".format(field, param_name))
            fig, ax = plt.subplots()
            norm = matplotlib.colors.Normalize(vmin=bound[0], vmax=bound[1])
            for item in data:
                param_val = item["parameters"][param_name]
                color = cmap(norm(param_val))
                handler(ax, item, color=color, lw=0.7)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            plt.colorbar(sm, label=param_name)

            ax.grid()
            dir_path = self.workspace.plots / "{}_by_param".format(field)
            dir_path.mkdir(exist_ok=True)
            path = dir_path / "{}.png".format(param_name)
            print("saving ({}, {}) to {}".format(field, param_name, path))
            fig.savefig(path, bbox_inches="tight", dpi=250)
            plt.close(fig)

    def _update_plots(self, cl_true, cl_nn, k_pk_true, pk_true, k_pk_nn, pk_nn):
        ell = cl_true["ell"]

        for _, ax in self.figs.values():
            ax.axhline(0, color="k")

        # TT
        tt_relerr = (cl_nn["tt"] - cl_true["tt"]) / cl_true["tt"]
        self.figs["tt"][1].semilogx(ell, tt_relerr, **LINESTYLE_RED)

        cosmic_variance = np.sqrt(2 / (2 * np.minimum(ell, 2000) + 1))
        tt_relerr_cv = tt_relerr / cosmic_variance
        self.figs["tt_cv"][1].semilogx(ell, tt_relerr_cv, **LINESTYLE_RED)

        def plot_err(ax, qty, style):
            ax.semilogx(ell, ell * (ell + 1) * qty, **style)

        def plot_err_pm(ax, qty, style):
            plot_err(ax, qty, style)
            plot_err(ax, -qty, style)

        # TE + EE
        for q in ("ee", "te"):
            err = (cl_nn[q] - cl_true[q])
            # err_relmax = err / cl_true[q].max()
            ax = self.figs[q][1]
            plot_err_pm(ax, cl_true[q] /  100.0, LINESTYLE_BLUE)
            plot_err_pm(ax, cl_true[q] / 1000.0, LINESTYLE_GREEN)
            plot_err(ax, err, LINESTYLE_RED)

        # P(k)
        # since P_NN(k) and P_true(k) may be sampled on different k grids, we
        # need to interpolate (in this case, onto the k_pk_true)
        pk_spline = CubicSpline(k_pk_nn, pk_nn)
        pk_nn_resampled = pk_spline(k_pk_true)
        pk_relerr = (pk_nn_resampled - pk_true) / pk_true
        self.figs["pk"][1].semilogx(k_pk_true, pk_relerr, **LINESTYLE_RED)
        self.figs["pk"][1].set_yscale("symlog", linthresh=0.01)

        # this will raise a warning that there are no positive values and
        # hence loglog is not possible in matplotlib version 3.3.1 but this
        # appears to be a bug in matplotlib and therefore can be ignored.
        self.figs["pk_abs"][1].loglog(k_pk_nn, pk_nn, **LINESTYLE_RED)
        self.figs["pk_abs"][1].loglog(k_pk_true, pk_true, **LINESTYLE_GREEN)

    def _save(self, prefix=None):
        for name, (fig, _) in self.figs.items():
            if not prefix:
                path = self.workspace.plots / "{}.png".format(name)
            else:
                dir_path = self.workspace.plots / prefix
                dir_path.mkdir(parents=True, exist_ok=True)
                path = dir_path / "{}.png".format(name)

            print("saving plot to", path)
            fig.savefig(path, dpi=200, bbox_inches="tight")

    def _close_figs(self):
        for (fig, _) in self.figs.values():
            plt.close(fig)
        self.figs = dict()
