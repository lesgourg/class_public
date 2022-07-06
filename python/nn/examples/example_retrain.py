"""
This script is written to manually parse arguments of the retrain (or fully train) process. The arguments are:
- w:    Path of the workspace where the networks and training data are to be stored. Per default it will point to ../../../classnet_workspace
- N:    List of networks to be retrained.
        select from ['Net_ST0_Reco','Net_ST0_Reio','Net_ST0_ISW','Net_ST1','Net_ST2_Reco','Net_ST2_Reio','Net_phi_plus_psi']
- g:    Provides a generation nametag. The nametag is expected to an interger. 
- s:    Steps of training process which are to be carried out. It expects an interger which consists of the didgets 1,2 and 3. These encode following steps:
        - 1:    Reading in networks. Read in covmat and bestfit to sample for the training, test and validationset. It is done in less than a minuit.
        - 2:    Generate sources from the sampled datasets. This takes several hours and should be parellized. Afterwards the k-array and the normalizations 
                for the networks are determined and stored within the function 'write_manifest'.
        - 3:    The training process is carried out. The training and model parameters are coded in the individual network classes under ~/path/to/class/python/nn/models/.
                It takes several hours and should be parallelized as well. Currently it is only available under CPU
- p:    Path to the covmat and bestfit file. These files should de denoted with {}.covmat and {}.bestfit. 
        If no path is provided the data folder in the workspace directory is used.

"""



from multiprocessing.sharedctypes import Value
import os
import sys
from classynet.workspace import Workspace, GenerationalWorkspace
import argparse
from pathlib import Path

import classy
import numpy as np

# Write argument parser
parser = argparse.ArgumentParser(
    prog='CLASSNet network retrain',
    description='This script is written to manually parse arguments of the retrain (or fully train) process'
)
parser.add_argument('-w', 
    help="Path of the workspace where the networks and training data are to be stored. Per default it will point to ../../../classnet_workspace", 
    default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../classnet_workspace"), 
    type=str)
parser.add_argument('-p', 
    help="Path to the covmat and bestfit file. These files should de denoted with {}.covmat and {}.bestfit. If no path is provided the data folder in the workspace directory is used.", 
    default=None, 
    type=str)
parser.add_argument('-n',
    help="List of networks to be retrained. Select from ['Net_ST0_Reco','Net_ST0_Reio','Net_ST0_ISW','Net_ST1','Net_ST2_Reco','Net_ST2_Reio','Net_phi_plus_psi']. Per default all networks are retrained.", 
    action='append', 
    nargs='+')
parser.add_argument('-g', 
    help="Provides a generation nametag. The nametag is expected to an interger. ",
    default=None,
    type=int)
parser.add_argument('-s', 
    help="Steps of training process which are to be carried out. 1: create/sample datasets. Computationally cheap. 2: run class full to obtain the source function. 3: Train the networks.",
    default=123,
    type=int)

args = parser.parse_args()

# write steps into list
steps = [int(d) for d in str(args.s)]

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

WORKSPACE_DIR = os.path.expanduser(args.w)

# Select a list of networks withgeneration name tag
full_network_list = ['Net_ST0_Reco','Net_ST0_Reio','Net_ST0_ISW','Net_ST1','Net_ST2_Reco','Net_ST2_Reio','Net_phi_plus_psi']

if args.n is None:
    my_network_list = full_network_list
else:
    my_network_list = args.n[0]
    for network in my_network_list:
        if network not in full_network_list:
            raise ValueError(network + " not known. Type -h for help.")

if args.g is not None:
    generations = {network: args.g for network in my_network_list}

    workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)
else:
    workspace = Workspace(WORKSPACE_DIR)

    for network in my_network_list:
        if os.path.exists(workspace.model_path(network)):
            raise ValueError("Trained network at " + str(workspace.model_path(network)) + " already exists. Either provide a nametag to train new networks or remove previous network from folder!")

assert isinstance(workspace, Workspace)

# Here we define the cosmological parameters on which the training is to be carried out
pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']

# Here we determine the used covmat and bestfit files.
if 1 in steps:
    if args.p is None:
        search_path = workspace.path / 'data'
    else:
        search_path = Path(os.path.expanduser(args.p))

    # list files in folder
    files = os.listdir(search_path)

    bf_file = [file for file in files if file[-8:]=='.bestfit']
    cm_file = [file for file in files if file[-7:]=='.covmat']

    if (len(bf_file)!=1) | (len(cm_file)!=1):
        raise ValueError("bestfit or covmat file not unambiguous")


    domain = workspace.domain_from_path(
        pnames         = pnames,
        bestfit_path   = search_path / bf_file[0],
        covmat_path    = search_path / cm_file[0],
        sigma_train    = 6,
        sigma_validation = 5,
        sigma_test     = 5,
    )

    #Save the domain within the workspace
    domain.save(workspace.domain_descriptor)

    #Sample parameters according to the domain and save them in the workspace for each dataset as "parameter_sample.h5"
    domain.sample_save(training_count=10000, validation_count=1000, test_count=1000)

if 2 in steps:
    #Load the data sets of parameters
    training, validation, test = workspace.loader().cosmological_parameters()
    
    # Generating training/validation/test data
    workspace.generator().generate_source_data(
        fixed=FIXED,
        training=training,
        validation=validation,
        test=test,
        fixed_training_only=FIXED_TRAINING_ONLY,
        processes=4)


    # Generate k array
    workspace.generator().generate_k_array()


    # Writing down the fixed and varying parameters which were used to generate the data
    workspace.generator().write_manifest(FIXED, training.keys())

if 3 in steps:
    # Import all networks
    from classynet.models import Net_ST0_Reco, Net_ST0_Reio, Net_ST0_ISW, Net_ST1, Net_ST2_Reco, Net_ST2_Reio, Net_phi_plus_psi

    NAME_FUNCTION_DICT = {'Net_ST0_Reco':Net_ST0_Reco, 'Net_ST0_Reio':Net_ST0_Reio, 'Net_ST0_ISW':Net_ST0_ISW, 'Net_ST1':Net_ST1, 'Net_ST2_Reco':Net_ST2_Reco, 'Net_ST2_Reio':Net_ST2_Reio, 'Net_phi_plus_psi':Net_phi_plus_psi}
    
    network_list = [NAME_FUNCTION_DICT[network] for network in my_network_list]
    
    # train all selected networks
    workspace.trainer().train_models(network_list, 8)

    # plot training progress
    workspace.plotter().plot_training_histories()
