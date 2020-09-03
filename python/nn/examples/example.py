import os
# from ..parameter_domain import EllipsoidDomain
from ..workspace import Workspace, GenerationalWorkspace
from ..parameter_sampling import DefaultParamDomain, EllipsoidDomain

# from classynet.parameter_domain import EllipsoidDomain
# from classynet.workspace import Workspace
import classy
import numpy as np

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

# # good generations, working with k sampling version 2
# generations = {
#     "Net_ST0_Reco":     10,
#     "Net_ST0_Reio":     6,
#     "Net_ST0_ISW":      6,
#     "Net_ST1":          6,
#     "Net_ST2_Reco":     6,
#     "Net_ST2_Reio":     6,
#     "Net_phi_plus_psi": 12
# }

# with mixed approach
# generations = {
#     "Net_ST0_Reco":     15,
#     "Net_ST0_Reio":     10,
#     "Net_ST0_ISW":      11,
#     "Net_ST1":          10,
#     "Net_ST2_Reco":     10,
#     "Net_ST2_Reio":     10,
#     "Net_phi_plus_psi": 17
# }


# generations = {
#     "Net_ST0_Reco":     18,
#     "Net_ST0_Reio":     13,
#     "Net_ST0_ISW":      14,
#     "Net_ST1":          13,
#     "Net_ST2_Reco":     13,
#     "Net_ST2_Reio":     13,
#     "Net_phi_plus_psi": 22,
# }

generations = {
    "Net_ST0_Reco":     20,
    "Net_ST0_Reio":     14,
    "Net_ST0_ISW":      15,
    "Net_ST1":          14,
    "Net_ST2_Reco":     14,
    "Net_ST2_Reio":     14,
    "Net_phi_plus_psi": 23,
}

workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)
assert isinstance(workspace, Workspace)

# domain = DefaultParamDomain(workspace.data / "base2018TTTEEE.covmat", sigma=5)

pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']
domain = EllipsoidDomain(
    workspace.data / "lcdm_11p_sn.bestfit",
    workspace.data / "lcdm_11p_sn.covmat",
    pnames=pnames,
    sigma=6,
)

# domain.sample_save(training_count=10000, validation_count=1000, path=workspace.data / "samples.h5")
# import sys; sys.exit(0)

training, validation = workspace.loader().cosmological_parameters()

# Generating training data
# workspace.generator().generate_data_for(
#     fixed=FIXED,
#     training=training,
#     validation=validation,
#     fixed_training_only=FIXED_TRAINING_ONLY,
#     processes=9)
# workspace.generator().generate_k_array()
# import sys; sys.exit(0)

# workspace.generator().write_manifest(FIXED, training.keys())

# workspace.generator().generate_data(FIXED DOMAIN, training=5, validation=2, processes=8)
# workspace.generator().generate_data(FIXED, DOMAIN, training=10000, validation=2000, processes=18)

# Training: any subset of models can be trained at once
# workspace.trainer().train_all_models(workers=36)

# from ..models import Net_ST0_Reco, Net_ST0_ISW
# workspace.trainer().train_models([
#     Net_ST0_Reco,
#     Net_ST0_ISW,
# ], 36)

# from ..models import Net_ST0_Reco
# workspace.trainer().train_model(Net_ST0_Reco, workers=18)
# import sys; sys.exit(0)

# from ..models import Net_phi_plus_psi
# workspace.trainer().train_model(Net_phi_plus_psi, workers=36)
# import sys; sys.exit(0)

# from ..models import Net_ST0_ISW
# workspace.trainer().train_model(Net_ST0_ISW, workers=18)
# import sys; sys.exit(0)

# ## Run CLASS for n cosmologies with and without NNs and produce error plots
tester = workspace.tester()
# tester.test(count=10000)

import matplotlib
matplotlib.use("agg")
plotter = workspace.plotter()
plotter.plot_spectra()
plotter.plot_source_functions()
plotter.plot_training_histories()
# plotter.plot_source_function_slice("t0_isw")
# triangle scatter plots of Cl/Pk errors vs. cosmological parameters
plotter.plot_scatter_errors()

# bm = workspace.benchmark_runner(warmup=1, iterations=50)
# bm.run(thread_counts=[1, 4])

# bm_plotter = workspace.benchmark_plotter()
# bm_plotter._load_data()
# bm_plotter.plot_perturbation_module()

# workspace.tester().test(50, processes=1, cheat=["t0_isw"], prefix="cheat_t0_isw")

# plotter.plot_training_histories()

# Compute Cl's with all source functions computed by CLASS _except_ one
if False:
    ALL_SOURCES = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m"]
    for i, select in enumerate(ALL_SOURCES):
        cheat = set(ALL_SOURCES) - set([select])
        workspace.tester().test(100, prefix=f"only_{select}", cheat=cheat)

