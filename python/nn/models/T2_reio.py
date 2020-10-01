import time
import os

import numpy as np
import h5py as h5

import torch
import torch.nn as nn
import torch.nn.functional as F

from .model import Model
from . import common
from .. import time_slicing

class Net_ST2_Reio(Model):

    HYPERPARAMETERS_DEFAULTS = {
            "learning_rate": 1e-3
            }

    def __init__(self, k, hp=None):
        super().__init__(k)

        n_inputs_cosmo = len(common.INPUTS_COSMO)
        n_inputs_tau = 1
        n_k = len(k)

        self.net_cosmo = nn.Linear(n_inputs_cosmo, 20)
        self.net_tau = nn.Linear(n_inputs_tau, 140)
        self.net_merge = nn.Sequential(
            nn.PReLU(),
            nn.Linear(20 + 140, 200),
            nn.PReLU(),
            nn.Linear(200, n_k),
        )

        if hp is None:
            hp = Net_ST2_Reio.HYPERPARAMETERS_DEFAULTS

        self.learning_rate = hp["learning_rate"]

    def forward(self, x):
        self.k_min = x["k_min"][0]
        y = self.net_merge(torch.cat((
            self.net_cosmo(common.get_inputs_cosmo(x)),
            self.net_tau(x["tau_relative_to_reio"][:, None]),
        ), axis=1))

        # linear_combination = self.net_basis(x)
        # correction = self.net_correction(x)
        # add = linear_combination + correction
        return x["g_reio"][:, None] * y

    def epochs(self):
        return 40

    def slicing(self):
        return time_slicing.TimeSlicingReio(0.6)

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

    def required_inputs(self):
        return set(common.INPUTS_COSMO + ["k_min", "tau_relative_to_reio", "g_reio"])

    def criterion(self):
        def loss(prediction, truth):
            return common.mse_truncate(self.k, self.k_min)(prediction, truth)
        return loss

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, lambda epoch: np.exp(-epoch / 8))

    def source_functions(self):
        return ["t2_reio"]
