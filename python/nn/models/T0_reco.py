import time
import os

import numpy as np
import h5py as h5

import torch
import torch.nn as nn
import torch.nn.functional as F

import pandas as pd

from .. import pytorch_spline

from .model import Model
# import classynet.models.common as common
from . import common
from .. import utils
from .. import time_slicing

PLOT_MODE = False

# TODO remove this in final
THESIS_MODE = False

def thesis_write(name, value):
    import os
    import pickle
    path = os.path.expanduser("~/masterarbeit/data/t0_reco.pickle")
    if os.path.exists(path):
        with open(path, "rb") as f:
            data = pickle.load(f)
    else:
        data = {}

    data[name] = value
    print("writing quantity", name)
    with open(path, "wb") as f:
        pickle.dump(data, f)


class BasisDecompositionNet(nn.Module):

    HYPERPARAMETERS_DEFAULTS = {
            "lin_cosmo_coeff": 20,
            "lin_tau_coeff": 80,
            "lin_coeff_1": 220,
            "lin_phases": 50,
            "relu_slope_basis": 0.3,
            }

    def __init__(self, k, n_inputs_cosmo, n_inputs_tau, n_k, hp=None):
        super().__init__()

        if hp is None:
            hp = BasisDecompositionNet.HYPERPARAMETERS_DEFAULTS

        self.k = k

        # 3x for cos, 2x for sin
        self.n_phases = 3 + 2
        # 1 damping envelope for cos, 1 for sin
        self.n_damping = 2
        # coefficients (cos, sin, psi_ref)
        self.n_coefficients = 2

        self.parameter_counts = [
                self.n_phases,
                self.n_damping,
                self.n_coefficients,
                ]
        self.parameter_count = sum(self.parameter_counts)

        self.n_spline_points = n_spline_points = 12
        self.spline_range = (self.k[0], 0.6)

        rtol = 1e-4
        self.k_spline = nn.Parameter(
                utils.powerspace((1 - rtol) * self.spline_range[0], (1 + rtol) * self.spline_range[1], 3, n_spline_points),
                # torch.linspace(self.spline_range[0], self.spline_range[1], n_spline_points),
                requires_grad=False
                )

        assert self.k_spline[0] <= self.spline_range[0]
        assert self.k_spline[-1] >= self.spline_range[1]

        self.spline_value_net_cosmo = nn.Linear(n_inputs_cosmo, 100)
        self.spline_value_net_tau = nn.Linear(1, 200)
        self.spline_value_net = nn.Sequential(
            nn.PReLU(),
            nn.Linear(100 + 200, 100),
            nn.PReLU(),
            nn.Linear(100, n_spline_points),
        )

        # TODO rename lin_cosmo_coeff because it is not only used for coefficients
        self.lin_cosmo = nn.Linear(n_inputs_cosmo, 50)
        self.lin_tau = nn.Linear(n_inputs_tau, 200)

        self.lin_parameters = nn.Sequential(
            nn.PReLU(),
            nn.Linear(self.lin_cosmo.out_features + self.lin_tau.out_features, 250),
            nn.PReLU(),
            nn.Linear(250, 100),
            nn.PReLU(),
            nn.Linear(100, self.parameter_count),
        )

    def spline_parameters(self):
        yield from self.spline_value_net_cosmo.parameters()
        yield from self.spline_value_net_tau.parameters()
        yield from self.spline_value_net.parameters()

    def forward(self, x):
        k_d = x["k_d"]
        r_s = x["r_s"]

        cosmo = common.get_inputs_cosmo(x)
        tau_g = common.get_inputs_tau_reco(x)

        # compute parameters
        parameter_inputs = torch.cat((
            self.lin_cosmo(cosmo),
            self.lin_tau(tau_g),
            ), dim=1)
        parameters = self.lin_parameters(parameter_inputs)

        phases, delta_k_d, coefficients = torch.split(parameters, self.parameter_counts, dim=1)

        # BASIS
        basis = self.basis(self.k, r_s, phases)

        k2 = self.k**2
        # DAMPING
        # k_over_k_d = self.k[None, :] / k_d[:, None]
        k_over_k_d2 = k2[None, :] / (k_d**2)[:, None]
        # k_over_k_d2 = k_over_k_d**2
        arg = -k_over_k_d2[None, ...] * (1 + delta_k_d.T[:, :, None])**2
        damping = torch.exp(arg)

        # SPLINE
        spline_values = self.spline_value_net(
                torch.cat((
                    self.spline_value_net_cosmo(cosmo),
                    self.spline_value_net_tau(x["tau_relative_to_reco"][:, None])
                    ), dim=1)
                )

        logk = torch.log(self.k_spline)
        spline = pytorch_spline.CubicSpline(logk, spline_values)
        B = torch.zeros((len(x["tau_relative_to_reco"]), len(self.k))).to(self.k.device)
        mask = (self.k >= self.spline_range[0]) & (self.k <= self.spline_range[1])
        k_eval = self.k[mask]
        B[:, mask] = spline(torch.log(k_eval))

        result = torch.cat((
            coefficients.T[..., None] * basis * damping,
            B[None, :, :],
            ), dim=0)

        if THESIS_MODE:
            thesis_write("k", self.k.detach().cpu().numpy())
            thesis_write("tau_rel", 10**x["tau_relative_to_reco"].detach().cpu().numpy())
            thesis_write("spline", B)
            thesis_write("trig", coefficients.T[..., None] * basis)
            thesis_write("damping", damping)
            thesis_write("trig + damping", coefficients.T[..., None] * basis * damping)
            thesis_write("basis", result.sum(axis=0))

        if PLOT_MODE:
            contributions = coefficients.T[..., None] * basis * damping
            import matplotlib
            matplotlib.use("agg")
            import matplotlib.pyplot as plt

            i_tau = 170 + 15
            print("@ tau/tau_rec ~", 10**(x["tau"][i_tau] - x["tau"][170]))

            # plt.figure()
            # plt.suptitle("damping")
            # plt.semilogx(k_standard.K_STANDARD, torch.exp(-(k[i_tau]/k_d[i_tau])**2), label="CLASS")
            # plt.semilogx(k_standard.K_STANDARD, damping[i_tau], label="NN")
            # plt.semilogx(k_standard.K_STANDARD, step[i_tau], label="step function")
            # plt.legend()
            # plt.grid()
            # plt.savefig(PLOT_DIR + "damping.png", dpi=200)
            # low = step * linear_combination_low
            # high = envelope * linear_combination_high

            plt.figure(figsize=(12, 8))
            plt.xlabel("$k [Mpc^{-1}]$")
            plt.title("$\\tau = {:.3f} \\tau_{{rec}}$".format(10**x["tau_relative_to_reco"][i_tau]))
            k = self.k
            plt.semilogx(k, result[0, i_tau], label=r"$\cos(\varphi_0^{cos} + k r_s (1 + \delta r_s^{cos}) + \varphi_2^{cos} k^2)$")
            plt.semilogx(k, result[1, i_tau], label=r"$\sin(\varphi_0^{sin} + k r_s (1 + \delta r_s^{sin}) + \varphi_2^{sin} k^2)/(k r_s)$")
            ln, = plt.semilogx(k, B[i_tau], label="spline")
            plt.scatter(self.k_spline, spline_values[i_tau], color=ln.get_color())


            plt.semilogx(k, result.sum(dim=0)[i_tau], ls="-.", color="r", label="Basis Network prediction")
            plt.semilogx(RESULT_K, RESULT[i_tau] / 0.08479736355112114, color="b", ls="-.", label="CLASS solution")

            plt.legend(loc="lower right")
            plt.grid()
            plt.savefig(PLOT_DIR + "basis.png", dpi=200)
            plt.xscale("linear")
            plt.savefig(PLOT_DIR + "basis_linear.png", dpi=200)

            plt.figure()
            from plotting.plot_source_function import plot_source_function
            plot_source_function(
                plt.gca(),
                self.k,
                10**x["tau_relative_to_reco"].cpu().numpy(),
                B.T.cpu().numpy(),
                levels=50
            )
            plt.savefig(PLOT_DIR + "offset_spline.png", dpi=200)

        return result

    def norm_sin(self, n):
        if n == -1:
            # if n == 0, the approximation is sin(k*r_s)/k, which, for small k, is r_s.
            # The growth of r_s is at some point counteracted by the damping term exp(-(k/k_D)^2).
            # The maximum of the approximation is at a value of order ~500.
            return 500
        else:
            # otherwise, the effective power of k is >= 0 and the envelope
            # has the same maximum as the cosine approximations (shifted by one power)
            return self.norm_cos(n)

    def norm_cos(self, n):
        # The envelope of the approximations is k^n * trig(k * r_s) * exp(-(k/k_D)^2)
        # (trig can be either sin or cos; we consider only n >= 0 here).
        # Analytically, one finds that this envelope is largest at initial times
        # (because it decreases monotonically because k_D decreases monotonically).
        # On the initial time slice, the maximum of the envelope is analytically found
        # to be (k_(D,ini))^n * (n/2)^(n/2) * exp(-n/2).
        # Typically, k_D has a value of ~0.3
        # at initial times (which in the considered cases is usually ~0.8 tau_rec).
        return 0.3**n * (n/2)**(n/2) * np.exp(-n/2)

    def basis_sin(self, arg, delta_r_s, delta_phi):
        args = arg[:, :, None] + self.k[None, :, None] * delta_r_s + delta_phi
        result = self.k_powers_sin * torch.sin(args) / self.sin_norms[None, None, :]
        return result

    def basis_cos(self, arg, delta_r_s, delta_phi):
        return [torch.cos(arg + self.k * delta_r_s[..., 0] + delta_phi[..., 0]) / self.norm_cos(0)]

    def basis(self, k, r_s, values):
        k_ = k[None, :]
        k_2 = self.k[None, :]**2
        r_s_ = r_s[:, None]
        arg = k_ * r_s_
        arg_inv = (1 / k_) * (1 / r_s_)
        values = values[:, None, :]
        arg_cos = values[..., 0] + arg * (1 + values[..., 1]) + k_2 * values[..., 2]
        arg_sin = arg * (1 + values[..., 3]) + k_2 * values[..., 4]
        out = torch.stack((
            torch.cos(arg_cos),
            torch.sin(arg_sin) * arg_inv
        ))

        return out

class CorrectionNet(nn.Module):

    HYPERPARAMETERS_DEFAULTS = {
            "lin_cosmo_corr": 50,
            "lin_tau_corr": 240,
            "lin_corr_1": 410,
            "lin_corr_2": 200,
            "relu_slope_corr": 0.3,
            }

    def __init__(self, n_inputs_cosmo, n_inputs_tau, n_k, hp=None):
        super().__init__()

        if hp is None:
            hp = CorrectionNet.HYPERPARAMETERS_DEFAULTS

        self.lin_cosmo_corr = nn.Linear(n_inputs_cosmo, hp["lin_cosmo_corr"])
        self.lin_tau_corr = nn.Linear(n_inputs_tau, hp["lin_tau_corr"])

        self.lin_corr_1 = nn.Linear(hp["lin_cosmo_corr"] + hp["lin_tau_corr"], hp["lin_corr_1"])
        self.lin_corr_2 = nn.Linear(hp["lin_corr_1"], hp["lin_corr_2"])
        self.lin_corr_3 = nn.Linear(hp["lin_corr_2"], n_k)


        self.layers = [
                self.lin_cosmo_corr,
                self.lin_tau_corr,
                self.lin_corr_1,
                self.lin_corr_2,
                self.lin_corr_3,
        ]

        # self.relu = nn.LeakyReLU(hp["relu_slope_corr"])
        self.relu = nn.PReLU()
        self.init_with_zero()

    def init_with_zero(self):
        with torch.no_grad():
            self.lin_corr_3.weight.zero_()
            self.lin_corr_3.bias.zero_()

    def forward(self, x):
        cosmo = common.get_inputs_cosmo(x)
        tau_g = common.get_inputs_tau_reco(x)

        y_cosmo = self.relu(self.lin_cosmo_corr(cosmo))
        y_tau = self.relu(self.lin_tau_corr(tau_g))

        y = torch.cat((y_cosmo, y_tau), axis=1)
        y = self.relu(self.lin_corr_1(y))
        y = self.relu(self.lin_corr_2(y))
        y = self.lin_corr_3(y)

        return y

class Net_ST0_Reco(Model):

    HYPERPARAMETERS_DEFAULTS = {
            "learning_rate": 1e-3
            }

    def __init__(self, k, hp=None):
        super().__init__(k)

        n_inputs_cosmo = len(common.INPUTS_COSMO)
        n_inputs_tau = 4
        n_k = len(self.k)

        self.net_basis = BasisDecompositionNet(self.k, n_inputs_cosmo, n_inputs_tau, n_k, hp=hp)
        self.net_correction = CorrectionNet(n_inputs_cosmo, n_inputs_tau, n_k, hp=hp)

        if hp is None:
            hp = Net_ST0_Reco.HYPERPARAMETERS_DEFAULTS

        self.learning_rate = hp["learning_rate"]

        weight = torch.ones_like(self.k)
        weight[self.k < 5e-3] *= 10
        # rescale so that total weight is the same
        self.loss_weight = nn.Parameter(weight / weight.sum() * len(k), requires_grad=False)

    def forward(self, x):
        self.k_min = x["k_min"][0]

        linear_combination = self.net_basis(x)
        correction = self.net_correction(x)

        if PLOT_MODE:
            import matplotlib
            matplotlib.use("agg")
            import matplotlib.pyplot as plt
            from plotting.plot_source_function import plot_source_function

            fig, axes = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))

            plot_source_function(
                axes[0, 0],
                self.k,
                10**x["tau_relative_to_reco"].cpu().numpy(),
                linear_combination.sum(dim=0).cpu().numpy().T,
                levels=50,
                title="linear combination"
            )
            plot_source_function(
                axes[0, 1],
                self.k,
                10**x["tau_relative_to_reco"].cpu().numpy(),
                correction.cpu().numpy().T,
                levels=50,
                title="correction"
            )
            plot_source_function(
                axes[1, 0],
                self.k,
                10**x["tau_relative_to_reco"].cpu().numpy(),
                (linear_combination.sum(dim=0) + correction).cpu().numpy().T,
                levels=50,
                title="total"
            )

            fig.savefig(PLOT_DIR + "parts.png", dpi=300)

        result = torch.cat((linear_combination, correction[None, :, :]), dim=0)
        result = result.sum(dim=0)

        if THESIS_MODE:
            thesis_write("k_min", self.k_min.item())
            thesis_write("linear_combination", linear_combination.cpu().detach().numpy())
            thesis_write("correction", correction.cpu().detach().numpy())
            thesis_write("result", result.cpu().detach().numpy())

        # import matplotlib
        # matplotlib.use("agg")
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.semilogx(self.k.cpu().detach(), result[153])
        # plt.savefig("/home/samaras/deleteme/t0_reco_slice.png", dpi=250)

        return result
        # return linear_combination + correction

    def epochs(self):
        return 40

    def optimizer(self):
        basis_params = list(self.net_basis.parameters())
        spline_params = list(self.net_basis.spline_parameters())
        basis_params = list(set(basis_params) - set(spline_params))

        return torch.optim.Adam([
            {"params": basis_params},
            {"params": spline_params},
            {"params": self.net_correction.parameters()}
            ], lr=self.learning_rate)

    def required_inputs(self):
        return set(common.INPUTS_COSMO + [
            "k_min",
            "tau_relative_to_reco", "e_kappa",
            "r_s", "k_d", "tau_relative_to_reco", "g_reco", "g_reco_prime"
        ])

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, [
            # basis net parameters
            lambda epoch: np.exp(-epoch / 8),
            # spline parameters
            lambda epoch: np.exp(-epoch / 8) if epoch < 5 else 0,
            # correction net parameter
            lambda epoch: 0 if epoch < 5 else np.exp(-(epoch - 5) / 8)
            ]
            )

    def source_functions(self):
        return ["t0_reco_no_isw"]

    def slicing(self):
        return time_slicing.TimeSlicingReco(4)

    def criterion(self):
        """Returns the loss function."""
        # TODO self.loss_weight?
        # def loss(prediction, truth):
        #     return common.mse_truncate(self.k, self.k_min)(prediction, truth)
        # return loss
        def loss(prediction, truth):
            mask = self.k >= self.k_min
            return torch.mean((prediction - truth)[:, mask]**2)
        return loss
