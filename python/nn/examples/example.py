import os
import sys
from classynet.workspace import Workspace, GenerationalWorkspace
from classynet.parameter_sampling import EllipsoidDomain

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

WORKSPACE_DIR = os.path.expanduser("~/software/class_versions/class_net/workspaces/test_training")

# IF THIS IS GONE THEN UNDO IS FINISHED

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

# domain = DefaultParamDomain(workspace.data / "base2018TTTEEE.covmat", sigma=5)

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
domain.sample_save(training_count=10, validation_count=5, test_count=5)
#Load the data sets of parameters
training, validation, test = workspace.loader().cosmological_parameters()

#print(training)
sys.exit()

# Generating training data
'''workspace.generator().generate_source_data(
    fixed=FIXED,
    training=training,
    validation=validation,
    test=test,
    fixed_training_only=FIXED_TRAINING_ONLY,
    processes=2)
'''
# Generate k array
#workspace.generator().generate_k_array()

# Writing down the fixed and varying parameters which were used to generate the data
#workspace.generator().write_manifest(FIXED, training.keys())

# Training: any subset of models can be trained at once
workspace.trainer().train_all_models(workers=36)

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

'''

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

#plotter.plot_training_histories()

ALL_SOURCES = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m", "delta_cb"]
interest=[
        #"t0_reco_no_isw", 
        #"t0_reio_no_isw", 
        #"t0_isw", 
        #"t1",
        #"t2_reco", 
        #"t2_reio", 
        #"phi_plus_psi", 
        #"delta_m", 
        "delta_cb",
        ]
ylim1 = {"tt":0.02,"ee":0.02,"te":1e-12,"pp":0.02}
ylim2 = {"tt":0.01,"ee":0.01,"te":5e-13,"pp":0.01}
ylim3 = {"tt":0.005,"ee":0.005,"te":1e-13,"pp":0.005}
# Compute Cl's with all source functions computed by CLASS _except_ one
#mode="only"
mode="except"
nonlinear="halofit"
#nonlinear="linear"
if True:
    import matplotlib
    matplotlib.use("agg")
    for i, select in enumerate(interest):
        subspace = workspace.sub("{}_{}".format(mode,select))
        if mode == "only":
            cheat = set(ALL_SOURCES) - set([select])
        elif mode == "except":
            cheat=set([select])
        else:
            raise ValueError("specify mode")
        # subspace.tester().test(1000, cheat=cheat, seed=1337)    
        subspace.tester().test(96, cheat=cheat, seed=1337,nonlinear=nonlinear)
        plotter = subspace.plotter()
        plotter.plot_spectra(include_params=False, suffix=nonlinear,ylim=ylim1)
        plotter.plot_spectra(include_params=False, suffix=nonlinear,ylim=ylim2)
        plotter.plot_spectra(include_params=False, suffix=nonlinear,ylim=ylim3)
        #plotter.plot_source_functions()
'''
