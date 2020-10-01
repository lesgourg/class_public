import os
from ..parameter_domain import EllipsoidDomain
from ..workspace import Workspace, GenerationalWorkspace
# from classynet.parameter_domain import EllipsoidDomain
# from classynet.workspace import Workspace
import classy
import numpy as np
from classy import Class

import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt

## Here, all the fixed parameters are specified
FIXED = {
    "output": "tCl,pCl,lCl,mPk",

    "N_ncdm": 1,
    # "N_ur": 2.0328,
    # "m_ncdm": 0.06,
    "deg_ncdm": 3,
    "Omega_Lambda": 0,

    "P_k_max_1/Mpc": 100.0,
    "l_max_scalars": 3000,

    "compute damping scale": "yes",
}

# these fixed parameters are for the generation
# of the training set only; they don't need to be
# set when networks are evaluated.
FIXED_TRAINING_ONLY = {
    # to avoid interpolation artifacts at the edges
    "k_min_tau0": 1e-4,
}


# PLANCK = {
#         "omega_b":   (0.02242, 0.00014),
#         "omega_cdm": (0.11933, 0.00091),
#         "H0":        (67.66,   0.42),
#         "tau_reio":  (0.0561,  0.0071)
#         }
# ## The domain on which the NNs should be trained is specified as a dict: name -> (min, max)
# DOMAIN = {p: (mu - 5 * sigma, mu + 5 * sigma) for p, (mu, sigma) in PLANCK.items()}

# WORKSPACE_DIR = "/scratch/work/samaras/delete_me/example/"
WORKSPACE_DIR = os.path.expanduser("~/CLASSnet_HPC/")

# those are the ones with Omega_k set to 0
# generations = {
#     "Net_ST0_Reco":     9,
#     "Net_ST0_Reio":     5,
#     "Net_ST0_ISW":      5,
#     "Net_ST1":          5,
#     "Net_ST2_Reco":     5,
#     "Net_ST2_Reio":     5,
#     "Net_phi_plus_psi": 10
# }

generations = {
    "Net_ST0_Reco":     10,
    "Net_ST0_Reio":      6,
    "Net_ST0_ISW":       6,
    "Net_ST1":           6,
    "Net_ST2_Reco":      6,
    "Net_ST2_Reio":      6,
    "Net_phi_plus_psi": 11
}

PARAMS_BESTFIT = {
    "omega_b":    0.022433,
    "omega_cdm":  0.12062,
    "h":          0.62585,
    "tau_reio":   0.054783,
    "w0_fld":     -0.28377,
    "wa_fld":     -2.457,
    "N_ur":       0.042083,
    "omega_ncdm": 0.0010267,
    "Omega_k":    -0.0044679
}

workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)


def get_sources(params):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute(level=["perturb"])

    result = cosmo.get_sources()
    return cosmo, result

shared = FIXED.copy()
shared.update(PARAMS_BESTFIT)

params = shared.copy()
params_nn = shared.copy()
params_nn.update({"neural network path": workspace})

print("running class (full)")
cosmo, (src, k, tau) = get_sources(params)
print("running class (nn)")
cosmo_nn, (src_nn, k_nn, tau_nn) = get_sources(params_nn)


tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
tau_eval = tau_rec * 0.95

k_min = cosmo.k_min()

idx_tau = np.argmin(np.abs(tau - tau_eval))
idx_tau_nn = np.argmin(np.abs(tau_nn - tau_eval))
print("idx_tau_nn", idx_tau_nn)

srcfun = "t0"

print("creating plot")
plt.figure(figsize=(6, 8))
plt.subplot(211)
plt.semilogx(k, src[srcfun][:, idx_tau], label="true")
plt.semilogx(k_nn, src_nn[srcfun][:, idx_tau], label="nn")
plt.axvline(k_min, label="k_min", color="k", ls="-.")
plt.grid()
plt.legend()


S = src[srcfun][:, idx_tau]
S_nn = src_nn[srcfun][:, idx_tau_nn]

import scipy.interpolate
S_nn_resampled = scipy.interpolate.CubicSpline(k_nn, S_nn)(k)

residual = (S_nn_resampled - S)

plt.subplot(212, sharex=plt.gca())
plt.semilogx(k, S / 1000, c="b", label="S_true / 1000")
plt.semilogx(k, residual, c="r", label="residual")
plt.legend()
plt.grid()


print("saving plot")
plt.savefig("/home/samaras/deleteme/{}.png".format(srcfun), dpi=200)
