import numpy as np

import matplotlib.pyplot as plt
from .plot_source_function import plot_source_function

class BasisPlotter:
    def __init__(self, size, with_coefficients=True, with_phases=True):
        self.with_coefficients = with_coefficients
        self.with_phases = with_phases

        self.size = size
        ncols = 1 + int(with_coefficients) + int(with_phases)

        if with_coefficients:
            self.col_envelope = 0
            self.col_src = 1
            self.col_phases = 2 if with_phases else None
        else:
            self.col_envelope = None
            self.col_src = 0
            self.col_phases = 1 if with_phases else None

        self.fig, self.axes = plt.subplots(ncols=ncols, nrows=size, sharey=True, figsize=(5 * ncols, 4 * size), squeeze=False)
        self.fig.subplots_adjust(left=0.1, right=1-0.1, wspace=0, top=1-0.05, bottom=0.05)

        self._create_labels()

        # LOL 
        if size % 2 == 0:
            self.axes[0][0].invert_yaxis()

    def _create_labels(self):
        if self.with_coefficients:
            self.axes[-1, self.col_envelope].set_xlabel("Envelope")
        if self.with_phases:
            self.axes[-1, self.col_phases].set_xlabel("$\\Delta \\varphi / \\pi$")

    def plot_coefficients(self, tau_rel, coefficients):
        for ax_coeff, coeff in zip(self.axes[:, self.col_envelope], coefficients.transpose()):
            ax_coeff.plot(coeff, tau_rel, "k-")
            ax_coeff.grid()
            ax_coeff.axvline(0, c="k", ls="--")

        return self

    def plot_basis(self, k, tau_rel, basis, labels=None):
        labels = labels if labels is not None else [None] * self.size

        for ax_src, basis_fun, label in zip(self.axes[:, self.col_src], basis.transpose(2, 1, 0), labels):
            plot_source_function(ax_src, k, tau_rel, basis_fun, 
                    show_ylabel=False, show_xlabel=False, 
                    levels=50,
                    title=label
                    )
        return self

    def plot_function(self, tau_rel, function, index, **kwargs):
        self.axes[index][self.col_src].plot(function, tau_rel, **kwargs)
        return self

    def plot_functions(self, tau_rel, functions, **kwargs):
        for i, f in enumerate(functions.T):
            self.plot_function(tau_rel, f, i, **kwargs)
        return self

    def plot_function_all(self, tau_rel, function, **kwargs):
        for i in range(self.size):
            self.plot_function(tau_rel, function, i, **kwargs)
        return self

    def plot_phases(self, tau_rel, phases, **kwargs):
        assert self.with_phases
        for ax, phase_func in zip(self.axes[:, self.col_phases], phases.T):
            ax.plot(phase_func / np.pi, tau_rel, **kwargs)
            ax.grid()
        return self

    def legend(self):
        for ax in self.axes[:, self.col_src]:
            ax.legend()
        return self

    @property
    def figure(self):
        return self.fig

if __name__ == "__main__":
    data = np.load("t0reco.npz")
    basis = data["basis"]
    coeff = data["coeff"]
    tau_rel = np.power(10, data["tau_rel"])
    r_s = data["r_s"]
    delta_r_s = data["delta_r_s"]
    k_d = data["k_d"]
    k_d_corr = data["corr_delta_k_d"][:, 0]
    delta_phi = data["delta_phi"]

    k_d_eff = k_d / k_d_corr

    print(basis.shape)
    print(coeff.shape)

    import k_standard
    k = k_standard.K_STANDARD

    damp = "e^{-(\\alpha k/k_d)^2}"
    labels_sin = ["$k^{{{}}} \\sin(k (r_s + \Delta r_s) + \\Delta \\varphi) {}$".format(i - 1, damp) for i in range(4)]
    labels_cos = ["$k^{{{}}} \\cos(k (r_s + \Delta r_s) + \\Delta \\varphi) {}$".format(i, damp) for i in range(4)]
    labels = labels_sin + labels_cos
    labels_psi = ["$$\\psi {}$".format(damp)]

    r_s_corr_sin = (r_s[:, None] + delta_r_s[:, 0, :])[..., :4]
    r_s_corr_cos = (r_s[:, None] + delta_r_s[:, 0, :])[..., 4:8]
    k_s = 2 * np.pi / r_s

    plotter = BasisPlotter(8, with_coefficients=True, with_phases=True)
    plotter.plot_basis(k, tau_rel, basis, labels=labels)
    plotter.plot_coefficients(tau_rel, coeff)

    plotter.plot_function_all(tau_rel, k_s, ls="-", c="k", label="$k_s$")
    plotter.plot_function_all(tau_rel, k_d, ls="-", c="r", label="$k_d$")
    plotter.plot_function_all(tau_rel, k_d_eff, ls="-.", c="r", label="$k_{d, eff.}$")

    plotter.plot_functions(tau_rel, 2 * np.pi / r_s_corr_sin, ls="-.", c="k", label="$k_{s, eff.}$")

    plotter.plot_phases(tau_rel, delta_phi[:, 0, :8], c="k", ls="-")

    plotter.legend()

    print("Saving...")
    plotter.figure.savefig("t0_reco_basis.png", dpi=250)

    ################################################################################
    # rest
    ################################################################################

    data = np.load("t0reco_parts.npz")
    lc = data["lc"]
    corr = data["corr"]

    fig, axes = plt.subplots(2, figsize=(5, 8))
    plot_source_function(axes[0], k_standard.K_STANDARD, tau_rel, lc.T, levels=50, title="linear comb.")
    plot_source_function(axes[1], k_standard.K_STANDARD, tau_rel, corr.T, levels=50, title="correction")

    fig.savefig("t0_reco_parts.png", dpi=200)
