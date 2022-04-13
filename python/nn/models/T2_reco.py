import time
import os

import numpy as np
import h5py as h5

import torch
import torch.nn as nn
import torch.nn.functional as F

from classynet.models import common
from classynet.models.model import Model
from classynet.tools import time_slicing

class Net_ST2_Reco(Model):

    HYPERPARAMETERS_DEFAULTS = {
            "learning_rate": 1e-3
            }

    def __init__(self, k, hp=None):
        super().__init__(k)

        n_inputs_cosmo = len(common.INPUTS_COSMO)
        n_inputs_tau = 1
        n_k = len(k)

        self.k = k
        self.net_cosmo = nn.Linear(n_inputs_cosmo, 20)
        self.net_tau = nn.Linear(n_inputs_tau, 140)
        self.net_combined = nn.Sequential(
            nn.PReLU(),
            nn.Linear(20 + 140, 200),
            nn.PReLU(),
            nn.Linear(200, 350),
            nn.PReLU(),
            nn.Linear(350, n_k),
        )

        if hp is None:
            hp = Net_ST2_Reco.HYPERPARAMETERS_DEFAULTS

        self.learning_rate = hp["learning_rate"]
        
        self.output_normalization = nn.Parameter(torch.ones(1), requires_grad=False)

    def forward(self, x):
        self.k_min = x["k_min"][0]
        inputs_cosmo = common.get_inputs_cosmo(x)
        inputs_tau = x["tau_relative_to_reco"][:, None]
        

        #TODO SG: Put the sqrt(6) somewhere else
        sqrt_6 = torch.sqrt(torch.tensor([6]))
        return x["g_reco"][:, None] * self.net_combined(
            torch.cat((
                self.net_cosmo(inputs_cosmo),
                self.net_tau(inputs_tau)
            ), dim=1))

        # linear_combination = self.net_basis(x)
        # correction = self.net_correction(x)
        # add = linear_combination + correction

        PLOT_MODE = False
        if PLOT_MODE:
            import matplotlib as mpl
            mpl.use("Qt4Agg")
            import matplotlib.pyplot as plt
            from plotting.plot_source_function import plot_source_function

            tau_rel = 10**x["tau_relative_to_reco"].cpu().numpy()
            g_reco = x["g_reco"][:, None].cpu().numpy()

            fig, axes = plt.subplots(nrows=2, ncols=2)
            plot_source_function(
                axes[0, 0],
                k_standard.K_STANDARD,
                tau_rel,
                (g_reco * linear_combination.cpu().numpy()).T,
                levels=50
            )
            plot_source_function(
                axes[0, 1],
                k_standard.K_STANDARD,
                tau_rel,
                (g_reco * correction.cpu().numpy()).T,
                levels=50
            )
            plot_source_function(
                axes[1, 0],
                k_standard.K_STANDARD,
                tau_rel,
                (g_reco * add.cpu().numpy()).T,
                levels=50
            )

            plt.show()

        # return x["g_reco"][:, None] * add
        return add

    def forward_reduced_mode(self, x, k_min_idx):
        self.k_min = x["k_min"][0]
        inputs_cosmo = common.get_inputs_cosmo(x)
        inputs_tau = x["tau_relative_to_reco"][:, None]

        #TODO SG: Put the sqrt(6) somewhere else
        y = x["g_reco"][:, None] * self.net_combined(
            torch.cat((
                self.net_cosmo(inputs_cosmo),
                self.net_tau(inputs_tau)
            ), dim=1)) 
        y = y * self.output_normalization#torch.tensor([0.0018118434977621952])

        #Do interpolation here 
        k_th = 1e-3
        i_k_th = np.argmin(np.abs(self.k - k_th))
        y[:,:i_k_th] = torch.outer(y[:,i_k_th] , (self.k[:i_k_th] / k_th)**2)

        return torch.flatten(y[:,k_min_idx:])



    def epochs(self):
        return 40

    def slicing(self):
        return time_slicing.TimeSlicingReco(4)

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        return torch.optim.Adam([
            # {"params": self.net_basis.parameters()},
            # {"params": self.net_correction.parameters()}
            ], lr=self.learning_rate)

    def criterion(self):
        def loss(prediction, truth):
            return common.mse_truncate(self.k, self.k_min)(prediction, truth)
            # import numpy as np
            # index_k_min = max(np.searchsorted(self.k.cpu().numpy(), self.k_min.cpu().item()) - 1, 0)
            # diff = (prediction - truth)[:, index_k_min:]
            # diff_rel = diff / truth[:, index_k_min:]
            # k_th = 1e-3
            # k = self.k[index_k_min:]
            # return torch.mean(diff[:, k > k_th]**2) + 1e-3 * torch.mean(diff_rel[:, k <= k_th]**2)
        return loss

    def required_inputs(self):
        return set(common.INPUTS_COSMO + ["k_min", "tau_relative_to_reco", "g_reco"])

    def tau_training(self):
        return None

        # with h5.File(os.path.join(os.path.expandvars("$CLASSNET_DATA"), "tau_t0_reco.h5"), "r") as f:
        #     tau_training = f["tau"][()]
        # return tau_training

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, lambda epoch: np.exp(-epoch / 8))
        # return torch.optim.lr_scheduler.LambdaLR(optimizer, [
        #     lambda epoch: np.exp(-epoch / 5),
        #     lambda epoch: 0 if epoch < 3 else np.exp(-(epoch - 3) / 5)
        #     ])

    def source_functions(self):
        return ["t2_reco"]

if __name__ == "__main__":
    iface = interface.TrainingInterface(Net_ST2_Reco)
    iface.run()
