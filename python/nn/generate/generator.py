import json
import sys
#from . import generate_data
from .generate_sources_files import generate_data, generate_parameters_and_data
from .generate_spectra import generate_spectral_data

class Generator:
    def __init__(self, workspace):
        self.workspace = workspace

    def generate_source_data(self, fixed=None, training=None, validation=None, test=None, processes=None, fixed_training_only=None):
        """
        this method works with precomputed `domain` data. It also writes a down the fixed Classy parameters and the parameters of te network
        """

        if fixed is not None:
            self.write_manifest(fixed, training.keys())

        # this must happen AFTER writing the manifest!
        if fixed_training_only is not None:
            fixed = fixed.copy()
            fixed.update(fixed_training_only)

        if training:
            generate_data(
                training,
                fixed,
                self.workspace.training_data, processes=processes
            )

        if validation:
            generate_data(
                validation,
                fixed,
                self.workspace.validation_data, processes=processes
            )

        if test:
            generate_data(
                test,
                fixed,
                self.workspace.test_data, processes=processes
            )

    def generate_spectra_data(self, training=None, validation=None, test=None, processes=None, fixed_nn_only=None, use_nn=False):
        """
        [SG]: TODO description
        """

        # load fixed CLASS parameters which were also used to generate the training samples
        fixed = self.workspace.loader().manifest()['fixed']

        if fixed_nn_only:
            fixed_nn_only['neural network path']=str(self.workspace.path)
        else:
            fixed_nn_only = {'neural network path':str(self.workspace.path)}

        if use_nn:
            path_suffix = 'NN_cls.h5'
        else:
            path_suffix = 'FULL_cls.h5'

        if training:
            generate_spectral_data(
                self.workspace,
                training,
                fixed,
                self.workspace.training_data / path_suffix, processes=processes, fixed_nn_only=fixed_nn_only, use_nn=use_nn,
            )

        if validation:
            generate_spectral_data(
                self.workspace,
                validation,
                fixed,
                self.workspace.validation_data / path_suffix, processes=processes, fixed_nn_only=fixed_nn_only, use_nn=use_nn,
            )

        if test:
            generate_spectral_data(
                self.workspace,
                test,
                fixed,
                self.workspace.test_data / path_suffix, processes=processes, fixed_nn_only=fixed_nn_only, use_nn=use_nn
            )

    def generate_k_array(self):
        import glob
        import numpy as np
        import h5py as h5
        #from tqdm import tqdm
        import matplotlib.pyplot as plt

        files = glob.glob(str(self.workspace.training_data / "sources_*.h5"))

        print("Creating standard k array...")
        ks = []
        #for fn in tqdm(files):
        i=0
        for fn in files:
            print("{}/{}".format(i,len(files)))
            with h5.File(fn, "r") as f:
                k = f["sampling/k"][()]
                ks.append(k)
            i+=1

        k_mins = np.array([k[0] for k in ks])
        k_maxs = np.array([k[-1] for k in ks])
        extra = np.geomspace(k_mins.min(), k_mins.max(), 50)

        k_longest = max(ks, key=len)
        k = np.unique(np.sort(np.concatenate((extra, k_longest, np.array([k_maxs.max()])))))

        print("created k sampling of {} points with k.min() = {}, k.max() = {}.".format(
            len(k), k[0], k[-1]))
        print("saving to {}.".format(self.workspace.k))
        np.save(self.workspace.data / "k.npy", k)

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

    def write_manifest(self, fixed, varying_names):
        # Save the manifest declaring the fixed and variable inputs
        # the data has been generated for
        with open(self.workspace.manifest, "w") as dest:
            json.dump(self.manifest(fixed, varying_names), dest)

    def manifest(self, fixed, varying_names):
        # SG: add cosmological fixed parameters into the varying_names list
        #necessary_NN_parameters = ['']
        #fixed_cosmological_names = []
        #print(fixed)
        #print(varying_names)
        #sys.exit()
        return {
            "fixed": fixed,
            "cosmological_parameters": list(varying_names)
        }

