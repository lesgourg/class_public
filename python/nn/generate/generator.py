import json
import sys
import os
from .generate_sources_files import generate_data
from .generate_spectra import generate_spectral_data

class Generator:
    def __init__(self, workspace):
        self.workspace = workspace

    def generate_source_data(self, fixed=None, training=None, validation=None, test=None, processes=None, fixed_training_only=None):
        """
        this method works with precomputed `domain` data. It also writes a down the fixed Classy parameters and the parameters of te network
        """

        # If there is already a manifest written, it means that we are simulating data for a retraining. 
        # It is important to check that we did not alter the fixed CLASS parameters!
        if os.path.exists(self.workspace.manifest):
            fixed_manifest = self.workspace.loader().manifest()['fixed']
            for key in fixed_manifest:
                if fixed_manifest[key] != fixed[key]:
                    raise AttributeError("Error when generating source data. Fixed CLASS parameter do not correspond to the fixed manifest parameter!")
        else:
            # If there is no manifest yet. We can write it now
            self.write_to_manifest('fixed',fixed)
            self.write_to_manifest('cosmological_parameters', list(training.keys()))


        if fixed_training_only is not None:
            fixed_training = fixed.copy()
            fixed_training.update(fixed_training_only)

        if training:
            generate_data(
                training,
                fixed_training,
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

        # If there is no k-array yet in the manifest, it is to be determined and stored in the manifest
        if 'k' not in self.workspace.loader().manifest():
            k = self.generate_k_array()
            self.write_to_manifest('k',list(k))

        # If there is no normalization yet in the manifest, it will be stored in the manifest
        if 'normalization' not in self.workspace.loader().manifest():
            with open(self.workspace.training_data / "min_max.json","r") as file:
                normalization = json.load(file)
            self.write_to_manifest('normalization',normalization)

    def generate_spectra_data(self, training=None, validation=None, test=None, processes=None, fixed_nn_only=None, use_nn=False):
        """
        This function takes the parameter sets and create Cl spectra and matter power spectra and stores them 
        to analyse the performance of the NN.
        """

        # load fixed CLASS parameters which were also used to generate the training samples
        fixed = self.workspace.loader().manifest()['fixed']

        if fixed_nn_only:
            fixed_nn_only['workspace_path']=str(self.workspace.path)
        else:
            fixed_nn_only = {'workspace_path':str(self.workspace.path)}

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
        """
        This function reads in the simulated sources from class which are within the training data set and 
        determines a unified k array, which is then handed over the the neural networks.
        """

        import glob
        import numpy as np
        import h5py as h5

        files = glob.glob(str(self.workspace.training_data / "sources_*.h5"))

        print("Creating standard k array...")
        ks = []
        for i,fn in enumerate(files):
            print("{}/{}".format(i,len(files)))
            with h5.File(fn, "r") as f:
                k = f["sampling/k"][()]
                ks.append(k)

        k_mins = np.array([k[0] for k in ks])
        k_maxs = np.array([k[-1] for k in ks])

        extra = np.geomspace(k_mins.min(), k_mins.max(), 50)

        k_longest = max(ks, key=len)
        k = np.unique(np.sort(np.concatenate((extra, k_longest, np.array([k_maxs.max()])))))

        print("created k sampling of {} points with k.min() = {}, k.max() = {}.".format(
            len(k), k[0], k[-1]))
        
        return k


    def write_to_manifest(self, key, value):
        """
        In the manifest the fixed and variabel cosmological inputs are stored. Once created during first training they are not to be alternated again.
        Furthermore the constant quantities which are used by all networks are stored here. These are namly:
        1. The input normalizations.
        2. The k-array.
        """
        if os.path.exists(self.workspace.manifest):
            my_manifest = self.workspace.loader().manifest()
        else:
            my_manifest = {}

        my_manifest[key] = value
        with open(self.workspace.manifest, "w") as dest:
            json.dump(my_manifest, dest)