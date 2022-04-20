import time
import os

import numpy as np
import h5py as h5

import torch
from torch import pow, log
import torch.nn as nn
import torch.nn.functional as F

from classynet.models.model import Model
from classynet.models import common
from classynet.tools.utils import Timer

# import constants
from classynet.tools.utils import C, THETA_CMB

#
# Approximation for the source functions S_cb, S_m, S_{phi_plus_psi}.
# Further outlined in appendix A1 of "CosmicNet II" paper
#
def TFacc(kk, 
    z_d, 
    Omega_matter, 
    Omega_baryon, 
    Omega_ncdm, 
    Omega_lambda, 
    D, 
    H, 
    hubble, 
    rs_drag, 
    keq, 
    a_eq, 
    redshift):


    kk = kk[None, :]
    D = D[:, None]
    H = H[:, None]
    redshift = redshift[:, None]

    degen_hdm = torch.tensor(3.).to(H.device)

    # derived quantities
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

    z_drag = z_d
    y_drag = z_equality/(1.0+z_drag)

    sound_horizon_fit = rs_drag

    p_c = 0.25*(5.0-torch.sqrt(1+24.0*f_cdm))
    p_cb = 0.25*(5.0-torch.sqrt(1+24.0*f_cb))

    omega_denom = H*C/hubble
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
    return tf_cb, tf_cbnu

def geomspace(start, stop, *args, **kwargs):
    return torch.logspace(np.log10(start), np.log10(stop), *args, **kwargs)

def split_tau(tau, tau_mid=10000):
    """
    Given an array `tau`, keep only those points above `tau_mid` and replace
    all points below `tau_mid` with a logarithmic sampling between `tau[0]`
    and `tau_mid` with the same logarithmic density as above `tau_mid`.
    """
    n_points_after = int((tau > tau_mid).sum())
    decades_after_mid = torch.log10(tau[-1] / tau_mid)
    points_per_decade = n_points_after / decades_after_mid
    decades_before_mid = torch.log10(tau_mid / tau[0])
    points_before_mid = int(points_per_decade * decades_before_mid)

    tau_ = torch.cat((
        geomspace(tau.min().item(), tau_mid, points_before_mid).to(tau.device),
        tau[tau >= tau_mid]
    ))
    return tau_

def subsample_tau(tau, tau_):
    """
    Given two arrays `tau` and `tau_`, return an array of indices (call it `indices`) into `tau`
    such that `tau[indices[i]]` is the closest point in `tau` to `tau[i]`.
    """
    # return np.unique(np.argmin(np.abs(tau_[:, None] - tau[None, :]), axis=1))
    return torch.argmin(torch.abs(tau_[:, None] - tau[None, :]), axis=1)

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
                nn.Linear(300, 3*n_k),
                )

        if hp is None:
            hp = Net_phi_plus_psi.HYPERPARAMETERS_DEFAULTS

        self.learning_rate = hp["learning_rate"]

        # contains values for current batch, set in `forward` and used for loss
        self.current = {}

        self.output_normalization = nn.Parameter(torch.ones(1), requires_grad=False)

    def forward(self, x):
        timer = Timer()
        timer.start("forward")

        ## to be used in loss function to subsample tau
        self.raw_tau = x["raw_tau"]

        self.k_min = x["k_min"][0]
        raw_tau = x["raw_tau"]
        a_eq = x["a_eq"][0]
        k_eq = x["k_eq"][0]
        z = x["z"] # tau_size array
        H = x["H"] # tau_size array
        D = x["D"] # tau_size array

        # Calculate constant factor alpha
        index_MD = torch.argmin(torch.abs(z - 50.0)) #gives the index where z is 50
        tau_md = raw_tau[index_MD]
        D_md = x["D"][index_MD]
        alpha = tau_md**2 / 7.8 / D_md

        # perform correction for N_ur != 0
        N_ur = x["raw_cosmos/N_ur"][0]
        F = (1 + 0.2271 * (3.046 + N_ur)) / (1 + 0.2271 * 3.046)

        # derive required quantities
        h = x["raw_cosmos/h"][0] / torch.sqrt(F)
        Omega_b = x["raw_cosmos/omega_b"][0] / h**2 / F
        Omega_ncdm = x["raw_cosmos/omega_ncdm"][0] / h**2 / F
        Omega_cdm = x["raw_cosmos/omega_cdm"][0] / h**2 / F
        Omega_k = x["raw_cosmos/Omega_k"][0]
        rs_drag = x["rs_drag"][0]
        z_d = x["z_d"][0]

        Omega_m = Omega_b + Omega_cdm + Omega_ncdm
        Omega_Lambda = 1.0 - Omega_b - Omega_cdm - Omega_k

        # Calculate approximation
        timer.start("tfacc")
        approx_cb, approx = TFacc(self.k, z_d, Omega_m, Omega_b, Omega_ncdm, Omega_Lambda, D, H, h, rs_drag, k_eq, a_eq, z)
        timer.end("tfacc")

        timer.start("approx normalize")
        approx /= approx[:, [0]]
        approx_cb /= approx_cb[:,[0]]
        timer.end("approx normalize")

        timer.start("create approx_stack")
        approx_stack = torch.empty((len(x["tau"]), len(self.k), 3), device=self.k.device)
        timer.end("create approx_stack")

        timer.start("copy phi+psi")
        approx_stack[:, :, 0] = approx
        timer.end("copy phi+psi")

        timer.start("approx_delta_m")
        approx_delta_m = -alpha.item() * approx * \
                (self.k2 + 3.*Omega_k*(h/C)**2)
        timer.end("approx_delta_m")

        timer.start("copy delta_m")
        approx_stack[:, :, 1] = approx_delta_m
        timer.end("copy delta_m")

        timer.start("approx_delta_cb")
        approx_delta_cb = -alpha.item() * approx_cb * \
                (self.k2 + 3.*Omega_k*(h/C)**2)
        approx_stack[:, :, 2] = approx_delta_cb
        timer.end("approx_delta_cb")

        inputs_cosmo = common.get_fields(x, self.cosmo_inputs())
        tau = x["tau"]
        inputs_tau = tau[:, None]

        timer.start("concat lin_cosmo, lin_tau")
        y = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau),
            ), dim=1)
        timer.end("concat lin_cosmo, lin_tau")

        timer.start("net_merge_corr")
        correction = self.net_merge_corr(y)
        timer.end("net_merge_corr")

        correction = correction.view(len(tau), -1, 3)

        timer.start("compute result")
        result = (1. + correction) * approx_stack
        timer.end("compute result")

        timer.end("forward")

        return result

    def forward_reduced_mode(self, x, k_min_idx):
        timer = Timer()
        timer.start("forward")

        ## to be used in loss function to subsample tau
        self.raw_tau = x["raw_tau"]

        self.k_min = x["k_min"][0]
        raw_tau = x["raw_tau"]
        a_eq = x["a_eq"][0]
        k_eq = x["k_eq"][0]
        z = x["z"] # tau_size array
        H = x["H"] # tau_size array
        D = x["D"] # tau_size array

        # Calculate constant factor alpha
        index_MD = torch.argmin(torch.abs(z - 50.0)) #gives the index where z is 50
        tau_md = raw_tau[index_MD]
        D_md = x["D"][index_MD]
        alpha = tau_md**2 / 7.8 / D_md

        # perform correction for N_ur != 0
        N_ur = x["raw_cosmos/N_ur"][0]
        F = (1 + 0.2271 * (3.046 + N_ur)) / (1 + 0.2271 * 3.046)

        # derive required quantities
        h = x["raw_cosmos/h"][0] / torch.sqrt(F)
        Omega_b = x["raw_cosmos/omega_b"][0] / h**2 / F
        Omega_ncdm = x["raw_cosmos/omega_ncdm"][0] / h**2 / F
        Omega_cdm = x["raw_cosmos/omega_cdm"][0] / h**2 / F
        Omega_k = x["raw_cosmos/Omega_k"][0]
        rs_drag = x["rs_drag"][0]
        z_d = x["z_d"][0]

        Omega_m = Omega_b + Omega_cdm + Omega_ncdm
        Omega_Lambda = 1.0 - Omega_b - Omega_cdm - Omega_k

        # Calculate approximation
        timer.start("tfacc")
        approx_cb, approx = TFacc(self.k, z_d, Omega_m, Omega_b, Omega_ncdm, Omega_Lambda, D, H, h, rs_drag, k_eq, a_eq, z)
        timer.end("tfacc")

        timer.start("approx normalize")
        approx /= approx[:, [0]]
        approx_cb /= approx_cb[:,[0]]
        timer.end("approx normalize")

        timer.start("create approx_stack")
        approx_stack = torch.empty((len(x["tau"]), len(self.k), 3), device=self.k.device)
        timer.end("create approx_stack")

        timer.start("copy phi+psi")
        approx_stack[:, :, 0] = approx * self.output_normalization
        timer.end("copy phi+psi")

        timer.start("approx_delta_m")
        approx_delta_m = -alpha.item() * approx * \
                (self.k2 + 3.*Omega_k*(h/C)**2)
        timer.end("approx_delta_m")

        timer.start("copy delta_m")
        approx_stack[:, :, 1] = approx_delta_m * D[:,None]
        timer.end("copy delta_m")

        timer.start("approx_delta_cb")
        approx_delta_cb = -alpha.item() * approx_cb * \
                (self.k2 + 3.*Omega_k*(h/C)**2)
        approx_stack[:, :, 2] = approx_delta_cb * D[:,None]
        timer.end("approx_delta_cb")

        inputs_cosmo = common.get_fields(x, self.cosmo_inputs())
        tau = x["tau"]
        inputs_tau = tau[:, None]

        timer.start("concat lin_cosmo, lin_tau")
        y = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau),
            ), dim=1)
        timer.end("concat lin_cosmo, lin_tau")

        timer.start("net_merge_corr")
        correction = self.net_merge_corr(y)
        timer.end("net_merge_corr")

        correction = correction.view(len(tau), -1, 3)

        timer.start("compute result")
        result = (1. + correction) * approx_stack
        timer.end("compute result")

        output = torch.flatten(result[:,k_min_idx:,:], start_dim=0,end_dim=1)

        timer.end("forward")
        return output

    def epochs(self):
        return 25

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

    def criterion(self):
        def loss(prediction, truth):
            # return common.mse_rel_()(prediction, truth)
            return common.mse_rel_truncate(self.k, self.k_min)(prediction, truth)
        return loss

        # TODO TODO TODO
        iterations = 0
        # weight = self.loss_weight[None, :, None]**2
        def loss(prediction, truth):
            nonlocal iterations
            print("using custom loss for phi+psi")

            # return torch.mean(((prediction - truth) / truth)**2)

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

            idx_k_min = torch.searchsorted(self.k, self.k_min)

            p_phipsi = prediction[:, :, 0]
            p_dm     = prediction[:, :, 1]
            p_dcm    = prediction[:, :, 2]

            # interpolation parameter
            x = (k_min - self.k[idx_k_min - 1]) / (self.k[idx_k_min] - self.k[idx_k_min - 1])
            # interpolate predicted value of delta_m at k_min between the two points
            # that surround k_min
            p_dm_log = torch.log(-p_dm)
            p_dm_log_interp = x * p_dm_log[:, idx_k_min - 1] + (1 - x) * p_dm_log[:, idx_k_min]
            p_dm_interp = -torch.exp(p_dm_log_interp)
            
            # do the same for delta_cb
            p_dcb_log = torch.log(-p_dcb)
            p_dcb_log_interp = x * p_dcb_log[:, idx_k_min -1] + (1 - x) * p_dcb_log[:, ixd_k_min]
            p_dcb_interp = -torch.exp(p_dcb_loginterp)

            # similarly for phi+psi (albeit without log)
            p_phipsi_interp = x * p_phipsi[:, idx_k_min - 1] + (1 - x) * p_phipsi[:, idx_k_min]

            prediction_eff = prediction[self.k >= self.k_min]
            prediction_eff = torch.insert(prediction, 0, torch.stack((p_phipsi_interp, p_dm_interp, p_dcb_interp), axis=1), 1)

            truth_eff = truth[self.k >= self.k_min]
            # in the training data, the value for k_min can be found at [idx_k_min - 1] because
            # the preprocessing does nearest neighbor extrapolation
            truth_eff = torch.insert(truth_eff, 0, truth[:, idx_k_min - 1], 1)


            # subsample tau since we are not interested in having a super high accuracy at
            # recombination and reionization
            loss_mask = np.unique(
                subsample_tau(self.raw_tau, split_tau(self.raw_tau, tau_mid=10000))
            )
            return torch.mean(((prediction_eff - truth_eff) / truth_eff)**2)[loss_mask]
            # weight[-1] = 50
            # return torch.mean(weight[:, None, :] * mask[None, :, None] * (prediction - truth)**2)
            # return torch.mean(weight[:, None, :] * mask[None, :, None] * ((prediction - truth) / truth)**2)
            # return torch.mean(((prediction - truth) / truth)[:, mask, :]**2)
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
            "raw_H", "raw_D", "raw_Omega_m",
            "tau", "k_eq",
            "rs_drag",
            "z_d", "H",
            ])

    def tau_training(self):
        return None

    def lr_scheduler(self, optimizer):
        return torch.optim.lr_scheduler.LambdaLR(optimizer, lambda epoch: np.exp(-epoch / 5))

    def source_functions(self):
        return ["phi_plus_psi", "delta_m", "delta_cb"]
