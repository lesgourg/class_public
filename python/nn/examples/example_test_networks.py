"""
This example script tests the performance of the network
"""

import os
import sys
from classynet.workspace import Workspace, GenerationalWorkspace
from classynet.parameter_sampling import EllipsoidDomain

import classy
import numpy as np

import torch

WORKSPACE_DIR = os.path.expanduser("~/software/class_versions/class_net/workspaces/test_11p_delete")

#create workspace instance
workspace = Workspace(WORKSPACE_DIR)

#create tester instance 
Tester = workspace.tester()

#parameters we want to be sampled from
pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']

#generate cls according to parameter sampling (if not done yet) 
# can take ~2h
Tester.create_cls(pnames,N=1000)

#load cls of FULL and NN and create comparisson plots
Tester.compare_cls_full_net(N_lines=100)
