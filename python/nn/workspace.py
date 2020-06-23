import os
import json
from pathlib import Path
from functools import wraps

import numpy as np

from . import training
from . import plotter
from .generate import generator
from .testing import tester

def create_dir(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        path = func(*args, **kwargs)
        path.mkdir(parents=True, exist_ok=True)
        return path
    return wrapper

class Workspace:
    """
    Represents the workspace (corresponding to a directory on disk) in which
    training/validation data, models, logs, plots, etc. will be stored.
    This class wraps around a directory and endows it with some utility methods.
    """

    def __init__(self, path):
        self.path = Path(path)
        self.path.mkdir(parents=True, exist_ok=True)

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

    @property
    @create_dir
    def data(self):
        """
        path to standard k array
        """
        return self.path / "data"

    @property
    @create_dir
    def plots(self):
        """
        path to directory where plots are saved
        """
        return self.path / "plots"

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

    def cosmological_parameters(self):
        return self.data / "samples.npz"

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
        return plotter.SourceFunctionPlotter(self)


class GenerationalWorkspace(Workspace):
    def __init__(self, path, generations):
        super().__init__(path)
        self.generations = generations

    @property
    @create_dir
    def plots(self):
        g = self.generations
        suffix = "_".join("{}_{}".format(k, g[k]) for k in sorted(self.generations))
        return self.path / ("plots_" + suffix)

    def model_path(self, name):
        return self.models / "{}_{}.pt".format(name, self.generations[name])

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
        path = self.workspace.cosmological_parameters()
        parameter_names = self.domain_descriptor()["fields"]
        def to_dict(dataset):
            return dict(zip(parameter_names, dataset.T))
        data = np.load(path)
        return to_dict(data["training"]), to_dict(data["validation"])

    def domain_descriptor(self):
        path = self.workspace.domain_descriptor()
        with open(path) as f:
            d = json.load(f)
        d["bestfit"] = np.array(d["bestfit"])
        d["covmat"] = np.array(d["covmat"])
        return d
