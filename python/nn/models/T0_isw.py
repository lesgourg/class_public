import numpy as np
import os

import torch
import torch.nn as nn
import torch.nn.functional as F

import h5py as h5

from classynet.models import common
from classynet.models.model import Model

class Net_ST0_ISW(Model):

    def __init__(self, k):
        super().__init__(k)

        n_inputs_cosmo = len(common.INPUTS_COSMO)
        n_inputs_tau = 2
        n_k = len(k)

        self.lin_cosmo = nn.Linear(n_inputs_cosmo, 100)
        self.lin_tau = nn.Linear(2, 5 * 50)

        self.lin_combined = nn.Sequential(
            nn.PReLU(),
            nn.Linear(self.lin_cosmo.out_features + self.lin_tau.out_features, 300),
            nn.PReLU(),
            nn.Linear(300, n_k),
            nn.PReLU(),
            nn.Linear(n_k, n_k)
        )

        # self.lin_cosmo = nn.Linear(n_inputs_cosmo, 20)
        # self.lin_tau = nn.Linear(2, 133)

        # self.lin_combined = nn.Sequential(
        #     nn.PReLU(),
        #     nn.Linear(self.lin_cosmo.out_features + self.lin_tau.out_features, 199),
        #     nn.PReLU(),
        #     nn.Linear(199, n_k)
        # )

        k_np = k.cpu().detach().numpy()
        density = np.gradient(np.cumsum(np.ones_like(k_np)), np.log(k_np))
        weight = 1. / density
        weight = weight / weight.sum() * n_k
        self.loss_weight = nn.Parameter(torch.from_numpy(weight).float(), requires_grad=False)

        self.output_normalization = nn.Parameter(torch.ones(1), requires_grad=False)

        # loss_weight = torch.ones_like(k)
        # loss_weight[k < 5e-3] *= 100
        # loss_weight = loss_weight / loss_weight.sum() * len(k)
        # self.loss_weight = nn.Parameter(loss_weight, requires_grad=False)


    def criterion(self):
        """Returns the loss function."""
        # TODO loss weight?
        def loss(prediction, truth):
            # TODO bit hacky
            return common.mse_truncate(self.k, self.k_min)(
                torch.sqrt(self.loss_weight[None, :]) * prediction,
                torch.sqrt(self.loss_weight[None, :]) * truth
            )
        return loss

        # def loss(prediction, truth):
        #     return torch.mean(self.loss_weight[None, :] * (prediction - truth)**2)
        # return loss

    def forward(self, x):
        self.k_min = x["k_min"][0]
        inputs_cosmo = common.get_inputs_cosmo(x)
        inputs_tau = torch.stack((x["tau"], x["D"]), dim=1)
        acc = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau)
            ), dim=1)
        prediction = x["e_kappa"][:, None] * self.lin_combined(
            torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau)
            ), dim=1)
        )
        return prediction


    def forward_reduced_mode(self, x, k_min_idx):
        self.k_min = x["k_min"][0]
        inputs_cosmo = common.get_inputs_cosmo(x)
        inputs_tau = torch.stack((x["tau"], x["D"]), dim=1)
        acc = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau)
            ), dim=1)
        prediction = x["e_kappa"][:, None] * self.lin_combined(
            torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau)
            ), dim=1)
        )
        return torch.flatten(prediction[:,k_min_idx:] * self.output_normalization)#torch.tensor([0.0011569450645745045]))



    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)

    def epochs(self):
        return 40

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, [
            lambda epoch: np.exp(-epoch / 8),
            ])

    def required_inputs(self):
        return set(common.INPUTS_COSMO + [
            "k_min",
            "tau",
            "D",
            # "tau_relative_to_reco",
            "e_kappa",
            ])

    def tau_training(self):
        return None
        # with h5.File(os.path.join(os.path.expandvars("$CLASSNET_DATA"), "tau_array.h5"), "r") as f:
        #     tau_training = f["tau"][()]
        # return tau_training

    def source_functions(self):
        return ["t0_isw"]

if __name__ == "__main__":
    iface = interface.TrainingInterface(Net_ST0_ISW)
    iface.run()
