"""
This example script outlines how to generate training data, carry out the process of training and analyzes the training process!
"""



import os
import sys
from classynet.workspace import Workspace, GenerationalWorkspace

import classy
import numpy as np

import torch
print("cuda_version:",torch.version.cuda)
print("available:",torch.cuda.is_available())
print("backends:", torch.backends.cudnn.enabled)
print("torch_version:",torch.__version__)

## Here, all the fixed parameters are specified
FIXED = {
    "output": "tCl,pCl,lCl,mPk",

    # "N_ur": 2.0328,
    # "m_ncdm": 0.06,
    "non linear": "halofit",
    "lensing": "yes",
    "N_ncdm": 1,
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
    # precision parameters
    "tol_background_integration":     1.e-12,
    "tol_perturb_integration":        1.e-6,
    "reionization_optical_depth_tol": 1.e-5,
}

# PLANCK = {
#         "omega_b":   (0.02242, 0.00014),
#         "omega_cdm": (0.11933, 0.00091),
#         "H0":        (67.66,   0.42),
#         "tau_reio":  (0.0561,  0.0071)
#         }
# ## The domain on which the NNs should be trained is specified as a dict: name -> (min, max)
# DOMAIN = {p: (mu - 5 * sigma, mu + 5 * sigma) for p, (mu, sigma) in PLANCK.items()}

WORKSPACE_DIR = os.path.expanduser("~/software/class_versions/class_net/workspaces/test_training")

# Select a generation name tag
generations = {
    "Net_ST0_Reco":     201,
    "Net_ST0_Reio":     201,
    "Net_ST0_ISW":      201,
    "Net_ST1":          201,
    "Net_ST2_Reco":     201,
    "Net_ST2_Reio":     201,
    "Net_phi_plus_psi": 201,
}

workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)

assert isinstance(workspace, Workspace)

pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']

domain = workspace.domain_from_path(
    pnames         = pnames,
    #bestfit_path   = "~/path/to/lcdm_11p_sn.bestfit",
    #covmat_path    = "~/path/to/lcdm_11p_sn.covmat",
    sigma_train    = 6,
    sigma_validation = 5,
    sigma_test     = 5,
)

#Save the domain within the workspace
domain.save(workspace.domain_descriptor)

#Sample parameters according to the domain and save them in the workspace as "samples.h5"
domain.sample_save(training_count=10000, validation_count=1000, test_count=1000)

#Load the data sets of parameters
training, validation, test = workspace.loader().cosmological_parameters()
 
# Generating training/validation/test data
workspace.generator().generate_source_data(
    fixed=FIXED,
    training=training,
    validation=validation,
    test=test,
    fixed_training_only=FIXED_TRAINING_ONLY,
    processes=2)


# Generate k array
workspace.generator().generate_k_array()


# Writing down the fixed and varying parameters which were used to generate the data
workspace.generator().write_manifest(FIXED, training.keys())


# Training: any subset of models can be trained at once
workspace.trainer().train_all_models(workers=36)

# Allso possible to train individual networks:
'''
from classynet.models import Net_ST0_Reco, Net_ST0_Reio, Net_ST0_ISW, Net_ST1, Net_ST2_Reco, Net_ST2_Reio, Net_phi_plus_psi
workspace.trainer().train_models([
    Net_phi_plus_psi,
    Net_ST0_Reco,
    Net_ST0_Reio,
    Net_ST0_ISW,
    Net_ST1,
    Net_ST2_Reco,
    Net_ST2_Reio,
], 8)
'''
