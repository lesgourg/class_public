"""
This script investigates on source functions created by CLASSFULL and CLASSNET.
It plots both of these contributions and investigates on the differences.
"""

import os
import sys
from classynet.workspace import Workspace, GenerationalWorkspace

import numpy as np

import torch

WORKSPACE_DIR = os.path.expanduser("../../../CLASSnet_Workspace")

#create workspace instance
workspace = Workspace(WORKSPACE_DIR)

#create tester instance 
Tester = workspace.tester()

#list of source functions we are interested in
#my_source_functions = ['phi_plus_psi','delta_cb','delta_m']
my_source_functions = ['phi_plus_psi','delta_cb','delta_m','t0_isw','t0_reco_no_isw','t0_reio_no_isw','t1','t2_reco','t2_reio']

# Plot source functions for the set of parameters at the center of the domain with CLASSFULL and CLASSNET
# The plots are stored within the workspace ~/path/to/workspace/results/source_functions/
Tester.plot_source_functions(my_source_functions,use_nn = True)
Tester.plot_source_functions(my_source_functions,use_nn = False)

# This routine generates the source functions for CLASSFULL and CLASSNET and compares the differences
# The plots are stored within the workspace ~/path/to/workspace/results/source_functions/
Tester.plot_source_functions_difference(my_source_functions)
