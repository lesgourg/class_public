import numpy as np
import os

import torch
import torch.nn as nn
import torch.nn.functional as F

import h5py as h5

from classynet.models.model import Model
from classynet.models import common
from classynet import utils
from classynet import time_slicing

class Net_ST0_Reio(Model):

    def __init__(self, k):
        super().__init__(k)

        n_inputs_cosmo = len(common.INPUTS_COSMO)
        n_inputs_tau = 4
        n_k = len(k)

        self.lin_cosmo = nn.Linear(n_inputs_cosmo, 20)
        self.lin_tau = nn.Linear(n_inputs_tau, 133)

        self.lin_combined = nn.Sequential(
            nn.PReLU(),
            nn.Linear(self.lin_cosmo.out_features + self.lin_tau.out_features, 500),
            nn.PReLU(),
            nn.Linear(500, n_k)
        )
    def forward(self, x):
        self.k_min = x["k_min"][0]

        inputs_cosmo = common.get_inputs_cosmo(x)
        inputs_tau = torch.stack([
            x["tau_relative_to_reio"],
            x["g_reio"],
            x["g_reio_prime"],
            x["e_kappa"],
        ], axis=1)

        prediction = self.lin_combined(
            torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau)
            ), dim=1)
        )
        return prediction


    def forward_reduced_mode(self, x, k_min_idx):
        self.k_min = x["k_min"][0]

        inputs_cosmo = common.get_inputs_cosmo(x)
        inputs_tau = torch.stack([
            x["tau_relative_to_reio"],
            x["g_reio"],
            x["g_reio_prime"],
            x["e_kappa"],
        ], axis=1)

        prediction = self.lin_combined(
            torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau)
            ), dim=1)
        )
        return prediction[:,k_min_idx:]


    def epochs(self):
        return 40

    def slicing(self):
        return time_slicing.TimeSlicingReio(0.6)

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)

    def lr_scheduler(self, opt):
        return torch.optim.lr_scheduler.LambdaLR(opt, lambda epoch: np.exp(-epoch / 8))

    def required_inputs(self):
        return set(common.INPUTS_COSMO + [
            "k_min",
            "tau_relative_to_reio",
            "e_kappa",
            "g_reio", "g_reio_prime",
            # "t0_reio_approx_1",
            ])

    def tau_training(self):
        with h5.File(os.path.join(os.path.expandvars("$CLASSNET_DATA"), "tau_t0_reio.h5"), "r") as f:
            tau_training = f["tau"][()]
        return tau_training

    def source_functions(self):
        return ["t0_reio_no_isw"]

    def criterion(self):
        def loss(prediction, truth):
            return common.mse_truncate(self.k, self.k_min)(prediction, truth)
        return loss

if __name__ == "__main__":
    iface = interface.TrainingInterface(Net_ST0_Reio)
    iface.run()
