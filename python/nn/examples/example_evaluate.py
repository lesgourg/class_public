import os
from classynet.workspace import Workspace
import classy
import numpy as np

FIXED = {
    "output": "tCl,pCl,lCl,mPk",

    "N_ncdm": 1,
    "N_ur": 2.0328,
    "m_ncdm": 0.06,

    # to avoid interpolation artifacts at the edges
    "P_k_max_1/Mpc": 10.0,
    # "k_min_tau0": 0.01,

    "compute damping scale": "yes",
}

params = {
    "omega_b":   0.02242,
    "omega_cdm": 0.11933,
    "H0":        67.66,
    "tau_reio":  0.0561,

    "neural network path": os.path.expanduser("~/path/to/workspace"),
}
params.update(FIXED)

cosmo = classy.Class()
cosmo.set(params)
cosmo.compute()
