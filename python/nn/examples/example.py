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

    # "N_ur": 2.0328,
    # "m_ncdm": 0.06,
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
    "tol_background_integration":     1.e-3,
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

# WORKSPACE_DIR = "/scratch/work/samaras/delete_me/example/"
WORKSPACE_DIR = os.path.expanduser("~/CLASSnet_HPC/")

generations = {
    "Net_ST0_Reco":     32,
    "Net_ST0_Reio":     23,
    "Net_ST0_ISW":      23,
    "Net_ST1":          22,
    "Net_ST2_Reco":     24,
    "Net_ST2_Reio":     23,
    "Net_phi_plus_psi": 31,
}

workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)

assert isinstance(workspace, Workspace)

# domain = DefaultParamDomain(workspace.data / "base2018TTTEEE.covmat", sigma=5)

pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']
domain = EllipsoidDomain.from_paths(
    bestfit_path   = workspace.data / "lcdm_11p_sn.bestfit",
    covmat_path    = workspace.data / "lcdm_11p_sn.covmat",
    pnames         = pnames,
    sigma_train    = 6,
    sigma_validate = 5,
)

# domain.save(workspace.domain_descriptor)
# domain.sample_save(training_count=10000, validation_count=1000, path=workspace.data / "samples.h5")
# import sys; sys.exit(0)

training, validation = workspace.loader().cosmological_parameters()

# # Generating training data
# workspace.generator().generate_data_for(
#     fixed=FIXED,
#     training=training,
#     validation=validation,
#     fixed_training_only=FIXED_TRAINING_ONLY,
#     processes=9)
# workspace.generator().generate_k_array()
# import sys; sys.exit(0)

# # Generating training data
# workspace.generator().generate_data_for(
#     fixed=FIXED,
#     training=None,
#     validation=validation,
#     fixed_training_only=FIXED_TRAINING_ONLY,
#     processes=9)
# import sys; sys.exit(0)

# workspace.generator().write_manifest(FIXED, training.keys())

# workspace.generator().generate_data(FIXED DOMAIN, training=5, validation=2, processes=8)
# workspace.generator().generate_data(FIXED, DOMAIN, training=10000, validation=2000, processes=18)

# Training: any subset of models can be trained at once
# workspace.trainer().train_all_models(workers=36)

# from ..models import Net_ST0_Reco, Net_ST0_ISW, Net_ST2_Reco, Net_ST2_Reio
# workspace.trainer().train_models([
    # Net_ST2_Reco,
    # Net_ST2_Reio,
# ], 36)

# from ..models import Net_ST0_Reco
# workspace.trainer().train_model(Net_ST0_Reco, workers=36)
# import sys; sys.exit(0)

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

ALL_SOURCES = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m"]

# ## Run CLASS for n cosmologies with and without NNs and produce error plots
workspace = workspace.sub("extrapolate t2_reco")
tester = workspace.tester()
tester.test(1000, seed=1234)
plotter = workspace.plotter()
plotter.plot_source_function_slice("t2_reco", marker=".", xlim=(1e-5, 1e-3))
plotter.plot_source_function_slice_tau("t2_reco")
plotter.plot_spectra(include_params=False)
# plotter.plot_source_function_slice("t2_reco")
# plotter.plot_source_functions()
# plotter.plot_scatter_errors()
# plotter.plot_training_histories()
import sys; sys.exit(0)

# bm = workspace.benchmark_runner(warmup=1, iterations=50)
# bm.run(thread_counts=[1, 4])

# bm_plotter = workspace.benchmark_plotter()
# bm_plotter._load_data()
# bm_plotter.plot_perturbation_module()

# workspace.tester().test(50, processes=1, cheat=["t0_isw"], prefix="cheat_t0_isw")

# plotter.plot_training_histories()

ALL_SOURCES = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m"]

# Compute Cl's with all source functions computed by CLASS _except_ one
if True:
    import matplotlib
    matplotlib.use("agg")
    for i, select in enumerate(["t0_isw"]):
        subspace = workspace.sub("only_{}".format(select))
        cheat = set(ALL_SOURCES) - set([select])
        subspace.tester().test(1000, cheat=cheat, seed=1337)
        plotter = subspace.plotter()
        plotter.plot_spectra(include_params=True)
        # plotter.plot_source_functions()
