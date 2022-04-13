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

# These parameters are to be adapted 
params = {
    "omega_b":   0.02238777,
    "omega_cdm": 0.1201992,
    "omega_ncdm": 0.0002748833,
    "h":        0.6891288,
    "tau_reio":  0.04941598,
    "Omega_Lambda": 0.0,
    "w0_fld": -0.9435546,
    "wa_fld": -0.2810646,
    "N_ur": 0.06651194,
    "Omega_k": -9.656155e-05,
    #"neural network path": os.path.expanduser("~/path/to/workspace"),
    "neural network path": os.path.expanduser("../../../CLASSnet_Workspace"),
    "use_nn": "true"
}

# merge both FIXED and adapted parameters
params.update(FIXED)

# call CLASS
cosmo = classy.Class()
cosmo.set(params)
cosmo.compute()

# read out spectra as usual using classy
my_cls = cosmo.lensed_cl(2000)

# create a plot if required
fig,ax = plt.subplots(figsize=(4,3))
ax.plot(my_cls['ell'][2:],my_cls['ell'][2:]*(my_cls['ell'][2:]+1)*my_cls['tt'][2:])
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell*(\ell+1) C_\ell^\mathrm{TT}$')
ax.grid(True)
plt.tight_layout()
fig.savefig('TT_spectrum.pdf')