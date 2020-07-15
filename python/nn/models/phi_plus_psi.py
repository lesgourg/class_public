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

def fitfunc(omega_matter, omega_baryon, omega_hdm, degen_hdm, omega_lambda, hubble, rs_drag, redshift, kk):
    """
/* This routine takes cosmological parameters and a redshift and sets up
all the internal scalar quantities needed to compute the transfer function. */
/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
                                in units of the critical density. */
/*        omega_baryon -- Density of baryons, in units of critical. */
/*        omega_hdm    -- Density of massive neutrinos, in units of critical */
/*        degen_hdm    -- (Int) Number of degenerate massive neutrino species */
/*        omega_lambda -- Cosmological constant */
/*        hubble       -- Hubble constant, in units of 100 km/s/Mpc */
/*        redshift     -- The redshift at which to evaluate */
/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
        sets many global variables for use in TFmdm_onek_mpc() */
        """
    from torch import log, sqrt

    def square(x):
        return x * x

    kk = kk[None, :]
    redshift = redshift[:, None]

    theta_cmb = 2.728/2.7; # Assuming T_cmb = 2.728 K

    # Look for strange input

    assert omega_baryon >= 0.0
    assert omega_hdm >= 0.0
    assert hubble > 0.0
    # assert torch.all(redshift > -1.0)
    # if redshift > 99.0:
    #       print("Large redshift entered.  TF may be inaccurate.")

    if degen_hdm < 1:
        degen_hdm = 1

    # Have to save this for TFmdm_onek_mpc()
    num_degen_hdm = degen_hdm;

    # This routine would crash if baryons or neutrinos were zero, so don't allow that
    if omega_baryon<=0:
        omega_baryon=1e-5;
    if omega_hdm<=0:
        omega_hdm=1e-5;

    omega_curv = 1.0-omega_matter-omega_lambda;
    omhh = omega_matter*square(hubble);
    obhh = omega_baryon*square(hubble);
    onhh = omega_hdm*square(hubble);
    f_baryon = omega_baryon/omega_matter;
    f_hdm = omega_hdm/omega_matter;
    f_cdm = 1.0-f_baryon-f_hdm;
    f_cb = f_cdm+f_baryon;
    f_bnu = f_baryon+f_hdm;

    # Compute the equality scale.
    z_equality = 25000.0*omhh/square(square(theta_cmb));        # Actually 1+z_eq
    k_equality = 0.0746*omhh/square(theta_cmb);

    # Compute the drag epoch and sound horizon
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*(1.0+z_drag_b1*pow(obhh,z_drag_b2));
    y_drag = z_equality/(1.0+z_drag);

    # G: changed to rs_drag
    # sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));
    sound_horizon_fit = rs_drag

    # Set up for the free-streaming & infall growth function
    p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

    omega_denom = omega_lambda+square(1.0+redshift)*(omega_curv+omega_matter*(1.0+redshift));
    omega_lambda_z = omega_lambda/omega_denom;
    omega_matter_z = omega_matter*square(1.0+redshift)*(1.0+redshift)/omega_denom;
    growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/ \
        (pow(omega_matter_z,4.0/7.0)-omega_lambda_z+ \
         (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
    growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0) \
                                                -omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
    growth_to_z0 = growth_k0/growth_to_z0;

    # Compute small-scale suppression
    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)* \
        pow(1+y_drag,p_cb-p_c)* \
        (1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/ \
        (1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))* \
        (1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = sqrt(alpha_nu);
    beta_c = 1/(1-0.949*f_bnu);
    # Done setting scalar variables
    hhubble = hubble;   # Need to pass Hubble constant to TFmdm_onek_hmpc()

    # ---------------------------- TFmdm_onek_mpc() ----------------------

    """
    Given a wavenumber in Mpc^-1, return the transfer function for the
    cosmology held in the global variables.
    Input: kk -- Wavenumber in Mpc^-1
    Output: The following are set as global variables:
            growth_cb -- the transfer function for density-weighted
                            CDM + Baryon perturbations.
            growth_cbnu -- the transfer function for density-weighted
                            CDM + Baryon + Massive Neutrino perturbations.
    The function returns growth_cb
    """
    qq = kk/omhh*square(theta_cmb);

    # Compute the scale-dependent growth functions
    y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))* \
                square(num_degen_hdm*qq/f_hdm);
    temp1 = pow(growth_k0, 1.0-p_cb);
    temp2 = pow(growth_k0/(1+y_freestream),0.7);
    growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
    growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

    # Compute the master function
    gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/(1+(kk*sound_horizon_fit*0.43)**4));
    qq_eff = qq*omhh/gamma_eff;

    tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
    tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
    tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*square(qq_eff));

    qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
    max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/ \
                (pow(qq_nu,-1.6)+pow(qq_nu,0.8));
    tf_master = tf_sup*max_fs_correction;

    # Now compute the CDM+HDM+baryon transfer functions
    tf_cb = tf_master*growth_cb/growth_k0;
    tf_cbnu = tf_master*growth_cbnu/growth_k0;
    # return tf_cb;
    # G: changed
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
        tau = x["tau"]
        k_eq = x["k_eq"]

        index_MD = torch.argmin(x["z"] - 10)
        tau_md = tau[index_MD]
        D_md = x["D"][index_MD]

        alpha = tau_md**2 / 10 / D_md

        N_ur = x["raw_cosmos/N_ur"][0]
        # perform correction for N_ur != 0
        F = (1 + 0.2271 * (3.046 + N_ur)) / (1 + 0.2271 * 3.046)

        h = x["raw_cosmos/h"][0] / torch.sqrt(F)
        Omega_b = x["raw_cosmos/omega_b"][0] / h**2 / F
        Omega_ncdm = x["raw_cosmos/omega_ncdm"][0] / h**2 / F
        Omega_cdm = x["raw_cosmos/omega_cdm"][0] / h**2 / F
        Omega_k = x["raw_cosmos/Omega_k"][0]
        rs_drag = x["rs_drag"][0]

        # approx = fitfunc_old(self.k, k_eq[0], Omega_b, Omega_cdm,h, rs_drag)
        Omega_m = Omega_b + Omega_cdm + Omega_ncdm
        Omega_Lambda = 1.0 - Omega_b - Omega_cdm - Omega_k
        approx = fitfunc(Omega_m, Omega_b, Omega_ncdm, 3, Omega_Lambda, h, rs_drag, x["z"], self.k)

        approx_k2 = approx * self.k2

        inputs_cosmo = common.get_fields(x, self.cosmo_inputs())
        inputs_tau = tau[:, None]

        y = torch.cat((
                self.lin_cosmo(inputs_cosmo),
                self.lin_tau(inputs_tau),
            ), dim=1)

        factors = self.net_merge_factor(y)

        # shape: (n_k, 2)
        approx_stack = torch.stack((approx, approx_k2), dim=2)

        correction = self.net_merge_corr(y)

        correction = correction.view(len(tau), -1, 2)

        # TODO CHANGE IF NORMALIZATION CHANGES!
        tau = self.current_tau = 10**tau

        # result_delta_m = (1 + correction[..., 1]) * factors[:, None, 1] * approx_k2
        # result_phi_plus_psi = (1 + correction[..., 0]) * factors[:, None, 0] * approx
        # result = torch.stack((result_phi_plus_psi, result_delta_m), dim=2)

        # correction = 1 + correction - torch.mean(correction, dim=1)[:, None, :]
        # result = correction * factors[:, None, :] * approx_stack[None, :, :]
        alpha_ = torch.Tensor([1.0, alpha.item()]).to(self.k.device)
        result = alpha_[None, None, :] * correction * approx_stack

        if False:
            # result_f_only = factors[:, None, :] * approx_stack[None, :, :]
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
            nonlocal iterations

            # if iterations % 1000 == 0:
            #     pdm0 = prediction[-1, :, 1]
            #     tdm0 = truth[-1, :, 1]
            #     import matplotlib; matplotlib.use("qt4agg")
            #     import matplotlib.pyplot as plt
            #     plt.semilogx(self.k.cpu().detach(), pdm0.cpu().detach(), label="pred")
            #     plt.semilogx(self.k.cpu().detach(), tdm0.cpu().detach(), label="truth")
            #     plt.axvline(self.k_min, c="k", ls="--", label="k_min")
            #     plt.grid()
            #     plt.legend()
            #     plt.show()
            # iterations += 1

            # return torch.mean((prediction - truth)**2)
            mask = self.k > self.k_min
            weight = torch.ones((truth.shape[0], 2)).to(truth.device)
            # weight[-1] = 50
            # return torch.mean(weight[:, None, :] * mask[None, :, None] * (prediction - truth)**2)
            return torch.mean(weight[:, None, :] * mask[None, :, None] * ((prediction - truth) / truth)**2)
        return loss

    def cosmo_inputs(self):
        return set(common.INPUTS_COSMO) - set(["cosmos/tau_reio"])

    def required_inputs(self):
        return self.cosmo_inputs() | set([
            "z",
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
