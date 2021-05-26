import os
# from ..parameter_domain import EllipsoidDomain
from classynet.workspace import Workspace, GenerationalWorkspace
from classynet.parameter_sampling import DefaultParamDomain, EllipsoidDomain

# from classynet.parameter_domain import EllipsoidDomain
# from classynet.workspace import Workspace
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
    
    #"input_verbose": 1,
    #"background_verbose":1,
    #"thermodynamics_verbose":1,
    #"perturbations_verbose":1,
    #"primordial_verbose":1,
    #"spectra_verbose":1,
    #"nonlinear_verbose":1,
    #"lensing_verbose":1,
    #"distortions_verbose":1,
    #"output_verbose":1,
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

#WORKSPACE_DIR = "/scratch/work/stadtmann/CLASSnet_Workspace/CLASSnet_Workspace_1/"
#WORKSPACE_DIR = os.path.expanduser("~/Class/CLASSnet_Workspace_old/")
#WORKSPACE_DIR = os.path.expanduser("~/Class/CLASSnet_Workspace_new/")
WORKSPACE_DIR = os.path.expanduser("~/Class/CLASSnet_Workspace_new_2/")

# IF THIS IS GONE THEN UNDO IS FINISHED

#generations = {
#    "Net_ST0_Reco":     32,
#    "Net_ST0_Reio":     23,
#    "Net_ST0_ISW":      23,
#    "Net_ST1":          22,
#    "Net_ST2_Reco":     24,
#    "Net_ST2_Reio":     23,
#    "Net_phi_plus_psi": 32,
#}
generations = {
    "Net_ST0_Reco":     101,
    "Net_ST0_Reio":     101,
    "Net_ST0_ISW":      101,
    "Net_ST1":          102,
    "Net_ST2_Reco":     101,
    "Net_ST2_Reio":     101,
    "Net_phi_plus_psi": 101,
}

workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)

assert isinstance(workspace, Workspace)

# domain = DefaultParamDomain(workspace.data / "base2018TTTEEE.covmat", sigma=5)

pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']

#domain = EllipsoidDomain.from_paths(
#    bestfit_path   = workspace.data / "lcdm_11p_sn.bestfit",
#    covmat_path    = workspace.data / "lcdm_11p_sn.covmat",
#    pnames         = pnames,
#    sigma_train    = 6,
#    sigma_validate = 5,
#)
#
#domain.save(workspace.domain_descriptor)
#domain.sample_save(training_count=1, validation_count=100, path=workspace.data / "samples.h5")

# import sys; sys.exit(0)

#training, validation = workspace.loader().cosmological_parameters()
# Generating training data

#workspace.generator().generate_data_for(
#    fixed=FIXED,
#    training=training,
#    validation=validation,
#    fixed_training_only=FIXED_TRAINING_ONLY,
#    processes=1)
#workspace.generator().generate_k_array()
# import sys; sys.exit(0)

# # Generating training data
# workspace.generator().generate_data_for(
#     fixed=FIXED,
#     training=None,
#     validation=validation,
#     fixed_training_only=FIXED_TRAINING_ONLY,
#     processes=9)
# import sys; sys.exit(0)

#workspace.generator().write_manifest(FIXED, training.keys())

# workspace.generator().generate_data_old(FIXED DOMAIN, training=5, validation=2, processes=8)
#workspace.generator().generate_data_old(FIXED, DOMAIN, training=10000, validation=1000, processes=32)

# Training: any subset of models can be trained at once
#workspace.trainer().train_all_models(workers=36)

#from classynet.models import Net_ST0_Reco, Net_ST0_Reio, Net_ST0_ISW, Net_ST1, Net_ST2_Reco, Net_ST2_Reio, Net_phi_plus_psi
#workspace.trainer().train_models([
    #Net_phi_plus_psi,
    #Net_ST0_Reco,
    #Net_ST0_Reio,
    #Net_ST0_ISW,
    #Net_ST1,
    #Net_ST2_Reco,
    #Net_ST2_Reio,
#], 8)

# from ..models import Net_ST0_Reco
# workspace.trainer().train_model(Net_ST0_Reco, workers=36)
#import sys; sys.exit(0)

# from ..models import Net_ST2_Reco
# workspace.trainer().train_model(Net_ST2_Reco, workers=36)

# from ..models import Net_ST0_ISW
# workspace.trainer().train_model(Net_ST0_ISW, workers=18)

# from ..models import Net_phi_plus_psi
# workspace.trainer().train_model(Net_phi_plus_psi, workers=36)
# import sys; sys.exit(0)

# from ..models import Net_ST0_ISW
# workspace.trainer().train_model(Net_ST0_ISW, workers=18)
# import sys; sys.exit(0)

#ALL_SOURCES = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m","delta_cb"]

# ## Run CLASS for n cosmologies with and without NNs and produce error plots
# workspace = workspace.sub("fix delta_m")
# tester = workspace.tester()
# tester.test(5, seed=1234)
#workspace = workspace.sub("pk_test")
plotter = workspace.plotter()
#plotter.plot_spectra(include_params=False)
# plotter.plot_source_function_slice("t2_reco", marker=".", xlim=(1e-5, 1e-3))
# plotter.plot_source_function_slice_tau("t2_reco")
# plotter.plot_source_function_slice("t2_reco")
# plotter.plot_source_functions()
# plotter.plot_scatter_errors()
# plotter.plot_training_histories()
#import sys; sys.exit(0)

#bm_plotter = workspace.benchmark_plotter()
#bm_plotter.plot_perturbation_module()
#import sys; sys.exit(0)

# workspace.tester().test(50, processes=1, cheat=["t0_isw"], prefix="cheat_t0_isw")

plotter.plot_training_histories()

#ALL_SOURCES = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m", "delta_cb"]
#interest=[
        #"t0_reco_no_isw", 
        #"t0_reio_no_isw", 
        #"t0_isw", 
        #"t1",
        #"t2_reco", 
        #"t2_reio", 
        #"phi_plus_psi", 
        #"delta_m", 
#        "delta_cb",
#        ]
# Compute Cl's with all source functions computed by CLASS _except_ one
#mode="only"
#mode="except"
#nonlinear="halofit"
#nonlinear="linear"
#if True:
#    import matplotlib
#    matplotlib.use("agg")
#    for i, select in enumerate(interest):
#        subspace = workspace.sub("{}_{}".format(mode,select))
#        if mode == "only":
#            cheat = set(ALL_SOURCES) - set([select])
#        elif mode == "except":
#            cheat=set([select])
#        else:
#            raise ValueError("specify mode")
#        # subspace.tester().test(1000, cheat=cheat, seed=1337)    
#        subspace.tester().test(96, cheat=cheat, seed=1337,nonlinear=nonlinear)
#        plotter = subspace.plotter()
#        plotter.plot_spectra(include_params=False, suffix=nonlinear)
#        #plotter.plot_source_functions()
