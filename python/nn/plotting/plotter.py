from classynet.plotting.source_function_plotter import SourceFunctionPlotter
from classynet.plotting.triangle_error_plotter import TrianglePlotter
from classynet.plotting.spectra_plotter import SpectraPlotter
from classynet.plotting.history_plotter import HistoryPlotter

class Plotter:

    def __init__(self, workspace):
        self.workspace = workspace

    def plot_training_histories(self):
        HistoryPlotter(self.workspace).plot_and_save()

    def plot_spectra(self, include_params=False,suffix=None,ylim=None):
        SpectraPlotter(self.workspace).plot(include_params=include_params,suffix=suffix,ylim=ylim)

    def plot_scatter_errors(self):
        TrianglePlotter(self.workspace).plot_and_save()

    def plot_source_functions(self):
        SourceFunctionPlotter(self.workspace).plot_source_functions()

    def plot_source_function_slice(self, *args, **kwargs):
        SourceFunctionPlotter(self.workspace).plot_slice(*args, **kwargs)

    def plot_source_function_slice_tau(self, *args, **kwargs):
        SourceFunctionPlotter(self.workspace).plot_slice_tau(*args, **kwargs)
