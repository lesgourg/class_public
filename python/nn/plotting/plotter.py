# from classynet.plotting.plot_source_function import SourceFunctionPlotter
# from classynet.plotting import TrianglePlotter
# from classynet.plotting import SpectraPlotter
# from classynet.plotting import HistoryPlotter

from classynet.plotting import plot_source_function

import numpy as np
import matplotlib.pyplot as plt

class Plotter:

    def __init__(self, workspace):
        self.workspace = workspace

    # def plot_training_histories(self):
    #     HistoryPlotter(self.workspace).plot_and_save()

    # def plot_spectra(self, include_params=False,suffix=None,ylim=None):
    #     SpectraPlotter(self.workspace).plot(include_params=include_params,suffix=suffix,ylim=ylim)

    # def plot_scatter_errors(self):
    #     TrianglePlotter(self.workspace).plot_and_save()

    # def plot_source_functions(self):
    #     SourceFunctionPlotter(self.workspace).plot_source_functions()

    # def plot_source_function_slice(self, *args, **kwargs):
    #     SourceFunctionPlotter(self.workspace).plot_slice(*args, **kwargs)

    # def plot_source_function_slice_tau(self, *args, **kwargs):
    #     SourceFunctionPlotter(self.workspace).plot_slice_tau(*args, **kwargs)

    def plot_k_array(self, k, k_mins):
        # create a histogram
        fig, ax = plt.subplots()
        bins = np.geomspace(k[0], k[-1], 50)
        ax.hist(k, bins=bins, histtype="step")
        ax.set_xscale("log")
        ax.set_xlabel("$k [Mpc^{-1}]$")
        fig.savefig(self.workspace.data / "k_hist.png", dpi=200)

        # create a scatter plot
        fig, ax = plt.subplots()
        ax.scatter(k, np.zeros_like(k), label="k_std")
        ax.scatter(k_mins, 1 + np.zeros_like(k_mins), label="k_mins")
        ax.set_xlim(1e-7, 2e-4)
        ax.set_xscale("log")
        ax.set_xlabel("$k [Mpc^{-1}]$")
        fig.savefig(self.workspace.data / "k_scatter.png", dpi=200)
