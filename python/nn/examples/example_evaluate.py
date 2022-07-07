import os
from classynet.workspace import Workspace
import classy
import numpy as np
import matplotlib.pyplot as plt

# This scipt is an example how to run CLASSNet for a set of parameters

# Fixed CLASS parameters. We require to ask for a massive neutrino species
FIXED = {
    "output": "tCl pCl lCl mPk",
    "N_ncdm": 1,
    "deg_ncdm": 3,
    'lensing': 'yes',
    "P_k_max_1/Mpc": 10.0,
    "compute damping scale": "yes",
}

# These parameters are to be adapted - now they are at the center of training domain
params = {
    "omega_b":   0.0224,
    "omega_cdm": 0.1216,
    "omega_ncdm": 0.000275,
    "h":        0.689,
    "tau_reio":  0.0494,
    "Omega_Lambda": 0.0,
    "w0_fld": -0.944,
    "wa_fld": -0.281,
    "N_ur": 0.066,
    "Omega_k": -9.66e-05,
    "A_s": 2.079e-9,
    "n_s": 0.971,
    "workspace_path": os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../classnet_workspace"),
    "use_nn": "yes",
    # decrease verbose if you prefer
    "nn_verbose":3
}

# merge both FIXED and adapted parameters
params.update(FIXED)

# call CLASS
cosmo = classy.Class()
cosmo.set(params)
cosmo.compute()

# read out spectra as usual using classy
my_cls = cosmo.lensed_cl(2000)
cosmo.struct_cleanup()

# create a plot if required
fig,ax = plt.subplots(figsize=(4,3))
ax.plot(my_cls['ell'][2:],my_cls['ell'][2:]*(my_cls['ell'][2:]+1)*my_cls['tt'][2:])
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell*(\ell+1) C_\ell^\mathrm{TT}$')
ax.grid(True)
plt.tight_layout()
fig.savefig('TT_spectrum.pdf')
