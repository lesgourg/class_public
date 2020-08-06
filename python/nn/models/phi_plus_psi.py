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

def fitfunc_old(k, keq, Omega_b, Omega0_cdm, h, rs_drag):
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



def TFacc(kk, z_d, Omega_matter, Omega_baryon, Omega_ncdm, Omega_lambda, D, H, hubble, rs_drag, keq, a_eq, redshift):
    theta_cmb = 2.7255/2.7

    kk = kk[None, :]
    D = D[:, None]
    H = H[:, None]
    redshift = redshift[:, None]

    # GS changed to 3
    degen_hdm = torch.tensor(3.).to(H.device)
    # degen_hdm = 1

    Omega_curv = 1.0-Omega_matter-Omega_lambda;
    omhh = Omega_matter*hubble**2
    obhh = Omega_baryon*hubble**2
    onhh = Omega_ncdm*hubble**2
    f_baryon = Omega_baryon/Omega_matter
    f_hdm = Omega_ncdm/Omega_matter
    f_cdm = 1.0-f_baryon-f_hdm
    f_cb = f_cdm+f_baryon
    f_bnu = f_baryon+f_hdm

    z_equality = 1./a_eq
    k_equality = keq

    from torch import pow, log
    z_drag = z_d
    y_drag = z_equality/(1.0+z_drag)

    sound_horizon_fit = rs_drag

    p_c = 0.25*(5.0-torch.sqrt(1+24.0*f_cdm))
    p_cb = 0.25*(5.0-torch.sqrt(1+24.0*f_cb))

    omega_denom = H*2997.92458/hubble
    growth_k0 = z_equality/D
    D0 = 1.0
    growth_to_z0 = z_equality/D0
    growth_to_z0 = growth_k0/growth_to_z0

    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*pow(1+y_drag,p_cb-p_c)*(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/(1-0.193*torch.sqrt(f_hdm*degen_hdm)+0.169*f_hdm*pow(degen_hdm,0.2))*(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = torch.sqrt(alpha_nu)
    gamma_eff = omhh*(alpha_gamma+(1-alpha_gamma)/(1+(0.43*rs_drag)**4))
    beta_c = 1/(1-0.949*f_bnu)

    qq = kk/keq/13.41


    y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*(degen_hdm*qq/f_hdm)**2
    temp1 = pow(growth_k0, 1.0-p_cb);
    temp2 = pow(growth_k0/(1+y_freestream),0.7);
    growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
    growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1

    gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/(1+(kk*sound_horizon_fit*0.43)**2))
    qq_eff = qq*omhh/gamma_eff;

    e = 2.718281828459045
    tf_sup_L = torch.log(e+1.84*beta_c*alpha_gamma*qq_eff);
    tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11))

    tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*(qq_eff)**2);

    qq_nu = 3.92*qq*torch.sqrt(degen_hdm/f_hdm);
    max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(degen_hdm,0.3+0.6*f_hdm)/(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
    tf_master = tf_sup*max_fs_correction;

    tf_cb = tf_master*growth_cb/growth_k0;
    tf_cbnu = tf_master*growth_cbnu/growth_k0;
    # return tf_cb, tf_cbnu
    return tf_cbnu

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
        raw_tau = x["raw_tau"]
        a_eq = x["a_eq"][0]
        k_eq = x["k_eq"][0]
        z = x["z"]
        H = x["H"]
        D = x["D"]

        index_MD = torch.argmin(torch.abs(z - 50.0))
        tau_md = raw_tau[index_MD]
        D_md = x["D"][index_MD]
        alpha = tau_md**2 / 7.8 / D_md

        N_ur = x["raw_cosmos/N_ur"][0]
        # perform correction for N_ur != 0
        F = (1 + 0.2271 * (3.046 + N_ur)) / (1 + 0.2271 * 3.046)

        h = x["raw_cosmos/h"][0] / torch.sqrt(F)
        Omega_b = x["raw_cosmos/omega_b"][0] / h**2 / F
        Omega_ncdm = x["raw_cosmos/omega_ncdm"][0] / h**2 / F
        Omega_cdm = x["raw_cosmos/omega_cdm"][0] / h**2 / F
        Omega_k = x["raw_cosmos/Omega_k"][0]
        rs_drag = x["rs_drag"][0]
        z_d = x["z_d"][0]

        # approx = fitfunc_old(self.k, k_eq[0], Omega_b, Omega_cdm,h, rs_drag)
        Omega_m = Omega_b + Omega_cdm + Omega_ncdm
        Omega_Lambda = 1.0 - Omega_b - Omega_cdm - Omega_k

        approx = TFacc(self.k, z_d, Omega_m, Omega_b, Omega_ncdm, Omega_Lambda, D, H, h, rs_drag, k_eq, a_eq, z)
        approx /= approx[:, [0]]

        # approx = fitfunc(Omega_m, Omega_b, Omega_ncdm, 3, Omega_Lambda, h, rs_drag, x["z"], self.k)

        # divide approximation by the SAME normalization constant as delta_m
        # TODO IMPORTANT WARNING DANGER do not hardcode this!
        approx_delta_m = -alpha.item() * approx / 154197.6664204682 * \
            (self.k2 + 3.*Omega_k*(h/2997.9)**2)

        inputs_cosmo = common.get_fields(x, self.cosmo_inputs())
        tau = x["tau"]
        inputs_tau = tau[:, None]

        y = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau),
            ), dim=1)

        factors = self.net_merge_factor(y)

        # shape: (n_k, 2)
        approx_stack = torch.stack((approx, approx_delta_m), dim=2)

        correction = self.net_merge_corr(y)

        correction = correction.view(len(tau), -1, 2)

        # TODO CHANGE IF NORMALIZATION CHANGES!
        tau = self.current_tau = 10**tau

        result = (1. + correction) * approx_stack

        if False:
            pdm0 = result[-1, :, 1].cpu().detach()
            k = self.k.cpu().detach()

            import matplotlib; matplotlib.use("qt4agg")
            import matplotlib.pyplot as plt

            plt.figure()
            plt.loglog(k, -pdm0, label="pred", marker="o")
            plt.axvline(self.k_min, c="k", ls="--", label="k_min")
            plt.grid()
            plt.legend()

            plt.show()


        return result

    def epochs(self):
        return 25

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

    def criterion(self):
        def loss(prediction, truth):
            return common.mse_rel_truncate(self.k, self.k_min)(prediction, truth)
        return loss

        # TODO TODO TODO
        iterations = 0
        # weight = self.loss_weight[None, :, None]**2
        def loss(prediction, truth):
            nonlocal iterations

            do_plot = False
            if do_plot and iterations % 1000 == 0:
                pdm0 = prediction[-1, :, 1].cpu().detach()
                tdm0 = truth[-1, :, 1].cpu().detach()
                k = self.k.cpu().detach()

                import matplotlib; matplotlib.use("qt4agg")
                import matplotlib.pyplot as plt

                plt.subplot(211)
                # plt.semilogx(self.k.cpu().detach(), pdm0.cpu().detach(), label="pred")
                # plt.semilogx(self.k.cpu().detach(), tdm0.cpu().detach(), label="truth")
                plt.loglog(k, -pdm0, label="pred")
                plt.loglog(k, -tdm0, label="truth")
                plt.axvline(self.k_min, c="k", ls="--", label="k_min")
                plt.grid()
                plt.legend()

                plt.subplot(212, sharex=plt.gca())
                rel_res = (pdm0 - tdm0) / tdm0
                plt.semilogx(k, rel_res)
                plt.grid()

                plt.show()

            iterations += 1

            return torch.mean(((prediction - truth) / truth)**2)
            # TODO the stuff below is not necessary, i think

            # Since in cases with Omega_k != 0, we don't go as low as
            # self.k[0], we need to treat those cases specially.
            # We cannot simply discard all points below self.k[0],
            # since during network evaluation, we are asked to return
            # a value for the source function at all k; including self.k[0].
            # Then, during evaluation, the value at k_min will be obtained by
            # means of interpolation.
            # This however means that we need to something similar in the
            # loss function in order to guarantee that those interpolated points
            # will be meaningful.
            # The strategy is the following:
            # If self.k[0] < k_min, we will discard points below k_min but
            # add an additional one at k_min by interpolation.
            # This way, the optimizer will tweak the network prediction
            # for k modes < k_min such that interpolation during evaluation
            # gives meaningful results.
            # The trick is the following: Since the training data is resampled
            # onto self.k using nearest neighbor extrapolation,
            # the value for self.k[0] will be incorrect but will be correct
            # for self.k_min.

            import ipdb; ipdb.set_trace()
            k_eff = self.k[self.k >= self.k_min]
            k_eff = torch.insert(k_eff, 0, self.k_min, 1)

            prediction_eff = prediction[self.k >= self.k_min]
            prediction_eff = torch.insert(prediction, 0, prediction[:, 0], 1)

            truth_eff = truth[self.k >= self.k_min]
            truth_eff = torch.insert(truth_eff, 0, truth[:, 0], 1)

            # return torch.mean((prediction - truth)**2)
            mask = self.k > self.k_min
            weight = torch.ones((truth.shape[0], 2)).to(truth.device)
            # weight[-1] = 50
            # return torch.mean(weight[:, None, :] * mask[None, :, None] * (prediction - truth)**2)
            # return torch.mean(weight[:, None, :] * mask[None, :, None] * ((prediction - truth) / truth)**2)
            return torch.mean(weight[:, None, :] * ((prediction - truth) / truth)[:, mask, :]**2)
        return loss

    def cosmo_inputs(self):
        inputs = common.INPUTS_COSMO[:]
        inputs.remove("cosmos/tau_reio")
        return inputs

    def required_inputs(self):
        return set(self.cosmo_inputs()) | set([
            "z",
            "k_min",
            "raw_tau",
            "a_eq",
            "raw_cosmos/h", "raw_cosmos/omega_b", "raw_cosmos/omega_cdm",
            "raw_cosmos/N_ur", "raw_cosmos/omega_ncdm", "raw_cosmos/Omega_k",
            "D",
            "tau", "k_eq",
            "rs_drag",
            "z_d", "H",
            ])

    def tau_training(self):
        return None

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, lambda epoch: np.exp(-epoch / 5))

    def source_functions(self):
        return ["phi_plus_psi", "delta_m"]
