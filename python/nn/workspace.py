import os
import json
from pathlib import Path
from functools import wraps

import numpy as np
import h5py as h5

from . import training
from . import plotter
from .generate import generator
from .testing import tester
from .benchmark import BenchmarkRunner
from .benchmark_plotter import BenchmarkPlotter

def create_dir(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        path = func(*args, **kwargs)
        path.mkdir(parents=True, exist_ok=True)
        return path
    return wrapper

class ResultDir:

    def __init__(self, root):
        self.root = root

    def sub(self, subdir):
        subdir = self.root / subdir
        subdir.mkdir(exist_ok=True)
        return ResultDir(subdir)

    @property
    @create_dir
    def plots(self):
        return self.root / "plots"

    @property
    def stats_file(self):
        return self.root / "errors.pickle"

class Workspace:
    """
    Represents the workspace (corresponding to a directory on disk) in which
    training/validation data, models, logs, plots, etc. will be stored.
    This class wraps around a directory and endows it with some utility methods.
    """

    def __init__(self, path, results=None):
        self.path = Path(path)
        self.path.mkdir(parents=True, exist_ok=True)

        if results is not None:
            self._results = results
        else:
            result_dir = self.path / "results"
            result_dir.mkdir(exist_ok=True)
            self._results = ResultDir(result_dir)

    def sub(self, sub):
        return Workspace(self.path, results=self._results.sub(sub))

    @property
    def plots(self):
        return self._results.plots

    @property
    def stats_file(self):
        return self._results.stats_file

    @property
    @create_dir
    def training_data(self):
        return self.path / "training"

    @property
    @create_dir
    def validation_data(self):
        return self.path / "validation"

    @property
    @create_dir
    def models(self):
        return self.path / "models"

    def model_path(self, name):
        return self.models / f"{name}.pt"

    def model_path_checkpoint(self, name, checkpoint):
        return self.models / f"{name}_checkpoint_{checkpoint}.pt"

    @property
    @create_dir
    def data(self):
        """
        path to standard k array
        """
        return self.path / "data"

    @property
    @create_dir
    def benchmark(self):
        """
        path to directory where benchmark results are saved
        """
        return self.plots / "benchmark"

    @property
    def benchmark_data(self):
        """
        path to data file inside benchmark directory
        """
        return self.benchmark / "data.json"

    @property
    def normalization_file(self):
        return self.training_data / "normalization.json"

    @property
    def manifest(self):
        return self.path / "manifest.json"

    @property
    def k(self):
        """
        path to standard k array
        """
        return self.data / "k.npy"

    @property
    @create_dir
    def history(self):
        return self.path / "history"

    def history_for(self, name):
        return self.history / f"{name}.csv"

    # def cosmological_parameters(self):
    #     return self.data / "samples.npz"

    def cosmological_parameters(self):
        return self.data / "samples.h5"

    def domain_descriptor(self):
        return self.data / "domain.json"

    def generator(self):
        return generator.Generator(self)

    def loader(self):
        return Loader(self)

    def trainer(self):
        return training.Trainer(self)

    def tester(self):
        return tester.Tester(self)

    def plotter(self):
        return plotter.Plotter(self)

    def benchmark_runner(self, warmup=2, iterations=50):
        return BenchmarkRunner(self, warmup=warmup, iterations=iterations)

    def benchmark_plotter(self):
        return BenchmarkPlotter(self)


class GenerationalWorkspace(Workspace):

    def __init__(self, path, generations, results=None):
        super().__init__(path, results=results)
        self.generations = generations

        if results is not None:
            self._results = results
        else:
            result_dir = self.path / "results"
            result_dir.mkdir(exist_ok=True)
            g = self.generations
            suffix = "_".join("{}_{}".format(k, g[k]) for k in sorted(self.generations))
            path = self.path / ("results_" + suffix)
            path.mkdir(exist_ok=True)
            self._results = ResultDir(path)

    def sub(self, sub):
        return GenerationalWorkspace(
            path=self.path,
            generations=self.generations,
            results=self._results.sub(sub))

    def model_path(self, name):
        return self.models / "{}_{}.pt".format(name, self.generations[name])

    def model_path_checkpoint(self, name, checkpoint):
        return self.models / "{}_{}_checkpoint_{}.pt".format(
            name, self.generations[name], checkpoint
        )

    def history_for(self, name):
        return self.history / f"{name}_{self.generations[name]}.csv"


class Loader:
    def __init__(self, workspace):
        self.workspace = workspace

    def manifest(self):
        with open(self.workspace.manifest) as src:
            return json.load(src)

    def k(self):
        return np.load(self.workspace.k)

    def cosmological_parameters(self):
        def load(g):
            return {key: value[()] for key, value in g.items()}
        with h5.File(self.workspace.cosmological_parameters(), "r") as f:
            return load(f["training"]), load(f["validation"])

    def domain_descriptor(self):
        path = self.workspace.domain_descriptor()
        with open(path) as f:
            d = json.load(f)
        d["bestfit"] = np.array(d["bestfit"])
        d["covmat"] = np.array(d["covmat"])
        return d

    def stats(self):
        """
        Load training stats.
        Requires that stats have already been generated by `Tester.test()`.
        """
        import pickle
        # TODO prefix?
        path = self.workspace.stats_file
        print("Loading stats from", path)
        with open(path, "rb") as f:
            return pickle.load(f)

    def benchmark_data(self):
        with open(self.workspace.benchmark_data) as f:
            return json.load(f)
