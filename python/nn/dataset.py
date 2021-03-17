import torch

import numpy as np
import h5py as h5
import os
import glob
from collections import defaultdict
import scipy
import scipy.interpolate
import functools
import multiprocessing
import random
from classy import Class

from classynet.data_providers import SourceFileDataProvider, CLASSDataProvider
from classynet import current_transformer

DEBUG_LOADING = False

class SourceFileDataset(torch.utils.data.Dataset):

    @staticmethod
    def for_model(model, directory, input_transformer=None, target_transformer=None, limit=None):
        if input_transformer is None:
            input_transformer = current_transformer.get_input_transformer_normalizer()
        if target_transformer is None:
            target_transformer = current_transformer.get_target_transformer_normalizer(k=k_standard.K_STANDARD)

        return SourceFileDataset(
                directory,
                input_transformer=input_transformer,
                target_transformer=target_transformer,
                source_function_selection=model.source_functions(),
                input_selection=model.required_inputs(),
                input_network_predictor_factory=model.input_networks if hasattr(model, "input_networks") else None,
                input_network_selection=model.input_network_selection() if hasattr(model, "input_network_selection") else None,
                k=k_standard.K_STANDARD,
                tau=model.tau_training(),
                limit=limit,
                )

    def __init__(self,
            data_directory,
            input_transformer, target_transformer,
            k=None,
            tau=None,
            source_function_selection=None,
            input_selection=None,
            limit=None):
        """
        :param data_directory: path of directory containing source files; source files must be named like "sources_{index}.h5".
        :type data_directory: int
        :param transform: Transformation instance
        :type transform: class:`transformation.Transformation`
        :param k: standard k array on which source functions are sampled
        :type k: class:`np.ndarray`
        :param tau: tau array on which to interpolate the source functions; defaults to specific tau array
                           given by CLASS for each cosmology
        :type tau: class:`np.ndarray`, optional
        :param source_function_selection: list of names of source functions to return, defaults to all source functions
        :type source_function_selection: [str], optional
        :param limit: maximum number of source files to load, useful for experimenting; defaults to full set of source files
        :type limit: int, optional
        :param batch_size: number of files per batch to return, defaults to 1
        :type batch_size: int, optional
        """

        self.k = k
        self.tau = tau

        file_names = glob.glob(os.path.join(data_directory, "sources_*.h5"))
        total_file_count = len(file_names)
        print("total file count in {}: {}".format(data_directory, total_file_count))
        # all_files = [os.path.join(data_directory, "sources_{}.h5".format(i)) for i in range(total_file_count)]
        all_files = file_names

        if limit:
            self.file_count = limit
            self.file_paths = random.sample(all_files, limit)
        else:
            self.file_count = total_file_count
            self.file_paths = all_files

        self.input_transformer = input_transformer
        self.target_transformer = target_transformer

        self.source_function_selection = source_function_selection
        self.input_selection = input_selection

    def __len__(self):
        return self.file_count

    def __getitem__(self, i):
        if DEBUG_LOADING:
            print("__getitem__({})".format(i))
        with h5.File(self.file_paths[i], "r") as source_file:
            return self._load_inputs_outputs(source_file)

    def _load_inputs_outputs(self, source_file):
        """
        Load, transform/normalize and return (inputs, outputs) for given source_file (h5py.File instance).
        """
        provider = SourceFileDataProvider(source_file, k_standard=self.k, compat_mode=False)

        raw_inputs = provider.get_inputs(tau=self.tau, k=self.k, input_selection=self.input_selection)
        raw_outputs = provider.get_outputs(tau=self.tau, k=self.k, source_function_selection=self.source_function_selection)

        inputs = self.input_transformer.transform_inputs(raw_inputs)
        outputs = self.target_transformer.transform_targets(raw_outputs, inputs=raw_inputs)

        inputs = self._to_tensor(inputs)
        outputs = self._to_tensor(outputs)

        if len(outputs) == 1:
            outputs = outputs[next(iter(outputs))]

        return inputs, outputs

    def _to_tensor(self, dct):
        return {k: torch.from_numpy(v).float() for k, v in dct.items()}



