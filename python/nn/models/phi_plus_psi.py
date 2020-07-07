import time
import os

import numpy as np
import h5py as h5

import torch
import torch.nn as nn
import torch.nn.functional as F

from .model import Model
from . import common
from .. import utils

def fitfunc(k, keq, Omega_b, Omega0_cdm, h, rs_drag):
  # Taken from http://background.uchicago.edu/~whu/transfer/tf_fit.c
    q = k/keq/13.41
    f_baryon = Omega_b / Omega0_cdm
    omhh = Omega0_cdm * h**2
    alpha_gamma = 1-0.328*torch.log(431.0*omhh)*f_baryon + 0.38*torch.log(22.3*omhh)*f_baryon**2
    gamma_eff = omhh*(alpha_gamma+(1-alpha_gamma)/(1+(0.43*rs_drag)**4))
    q_eff = q*omhh/gamma_eff
    q = q_eff
    T_0_L0 = torch.log(2.0*np.e+1.8*q)
    T_0_C0 = 14.2 + 731.0/(1+62.5*q)
    # return q*q*T_0_L0/(T_0_L0+T_0_C0*q*q)
    # corresponds to phi,psi and delta_m/k^2
    return T_0_L0/(T_0_L0+T_0_C0*q*q)

class Net_phi_plus_psi(Model):

    HYPERPARAMETERS_DEFAULTS = {
            "learning_rate": 1e-3
            }

    def __init__(self, k, hp=None):
        super().__init__(k)

        # kappa_reio not needed
        n_inputs_cosmo = len(common.INPUTS_COSMO) - 1
        n_k = len(k)

        self.logk = nn.Parameter(torch.log(self.k), requires_grad=False)
        self.k2 = self.k**2

        self.lin_cosmo = nn.Linear(n_inputs_cosmo, 20)
        self.lin_tau = nn.Linear(1, 100)

        self.net_merge_corr = nn.Sequential(
                nn.PReLU(),
                nn.Linear(20 + 100, 300),
                nn.PReLU(),
                nn.Linear(300, 2 * n_k),
                )

        self.net_merge_factor = nn.Sequential(
                nn.PReLU(),
                nn.Linear(20 + 100, 2),
                )

        if hp is None:
            hp = Net_phi_plus_psi.HYPERPARAMETERS_DEFAULTS

        self.learning_rate = hp["learning_rate"]

        # contains values for current batch, set in `forward` and used for loss
        self.current = {}

    def forward(self, x):
        self.k_min = x["k_min"][0]

        tau = x["tau"]
        k_eq = x["k_eq"]

        h = x["raw_cosmos/h"][0]
        Omega_b = x["raw_cosmos/omega_b"][0] / h**2
        Omega_cdm = x["raw_cosmos/omega_cdm"][0] / h**2
        rs_drag = x["rs_drag"][0]

        approx = fitfunc(self.k, k_eq[0], Omega_b, Omega_cdm,h, rs_drag)

        approx_k2 = approx * self.k2

        inputs_cosmo = common.get_fields(x, self.cosmo_inputs())
        inputs_tau = tau[:, None]

        y = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau),
            ), dim=1)

        factors = self.net_merge_factor(y)

        # shape: (n_k, 2)
        approx_stack = torch.stack((approx, approx_k2), dim=1)

        correction = self.net_merge_corr(y)

        correction = correction.view(len(tau), -1, 2)

        # TODO CHANGE IF NORMALIZATION CHANGES!
        tau = self.current_tau = 10**tau

        result = (1 + correction) * factors[:, None, :] * approx_stack[None, :, :]
        result_f_only = factors[:, None, :] * approx_stack[None, :, :]

        if False:
            import matplotlib
            matplotlib.use("qt5agg")
            import matplotlib.pyplot as plt
            plt.subplot(211)

            dm_class = -RESULT[-1, :, 1] / 22008.821725
            dm_net = result[-1, :, 1]
            from scipy.interpolate import CubicSpline
            dm_net_on_class_k = CubicSpline(self.k.cpu().detach(), dm_net)(RESULT_K)

            plt.loglog(RESULT_K, dm_class, label="$\\delta_m$ CLASS")
            plt.loglog(self.k.cpu().detach(), dm_net, label="$\\delta_m$ prediction: $(1+corr(k, \\tau)) \\alpha(\\tau) * fitfunc$")
            plt.loglog(self.k.cpu().detach(), result_f_only[-1, :, 1], label="$\\alpha(\\tau) * fitfunc$", ls="-.")
            plt.grid()

            plt.legend()

            plt.subplot(212, sharex=plt.gca())
            plt.ylabel("relative error")
            plt.axhline(0, color="C0", ls="--")
            plt.semilogx(RESULT_K, (dm_net_on_class_k - dm_class) / dm_class, label="relative error", color="C1")
            plt.xlabel("$k$")
            plt.legend()
            plt.grid()

            plt.show()


        return result
        # return basis + correction

    def epochs(self):
        return 20

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

    def criterion(self):
        iterations = 0
        # weight = self.loss_weight[None, :, None]**2
        def loss(prediction, truth):

            # nonlocal iterations
            # if iterations % 100 == 0:
            #     import matplotlib; matplotlib.use("qt4agg")
            #     import matplotlib.pyplot as plt
            #     plt.loglog(self.k.cpu().detach(), -prediction[-1, :, 1].cpu().detach(), label="pred")
            #     plt.loglog(self.k.cpu().detach(), -truth[-1, :, 1].cpu().detach(), label="truth")
            #     plt.axvline(self.k_min, c="k", ls="--", label="k_min")
            #     plt.grid()
            #     plt.legend()
            #     plt.show()
            # iterations += 1

            # return torch.mean((prediction - truth)**2)
            mask = self.k > self.k_min
            return torch.mean(mask[None, :, None] * ((prediction - truth) / truth)**2)
        return loss

    def cosmo_inputs(self):
        return set(common.INPUTS_COSMO) - set(["cosmos/tau_reio"])

    def required_inputs(self):
        return self.cosmo_inputs() | set([
            "k_min",
            "raw_cosmos/h", "raw_cosmos/omega_b", "raw_cosmos/omega_cdm",
            "D",
            "tau", "k_eq",
            "rs_drag",
            ])

    def tau_training(self):
        return None

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, lambda epoch: np.exp(-epoch / 5))

    def source_functions(self):
        return ["phi_plus_psi", "delta_m"]
