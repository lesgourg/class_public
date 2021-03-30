import json
from . import generate_data
from .generate_sources_files import generate_data, generate_parameters_and_data

class Generator:
    def __init__(self, workspace):
        self.workspace = workspace

    def generate_data(self, fixed, domain, training, validation, processes=None):
        self.write_manifest(fixed, domain.keys())

        generate_parameters_and_data(
            training,
            domain, fixed,
            self.workspace.training_data, processes=processes
        )
        generate_parameters_and_data(
            validation,
            domain, fixed,
            self.workspace.validation_data, processes=processes
        )

    def generate_data_for(self, fixed, training=None, validation=None, processes=None, fixed_training_only=None):
        """
        this method works with precomputed `domain` data.
        """

        if training:
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
        return {
            "fixed": fixed,
            "cosmological_parameters": list(varying_names)
        }
