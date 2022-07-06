"""
This example script tests the performance of the network by calculating the power spectra using CLASSNet and default CLASS
"""

import os
from classynet.workspace import Workspace, GenerationalWorkspace

WORKSPACE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../classnet_workspace")

#create workspace instance
workspace = Workspace(WORKSPACE_DIR)

#create tester instance 
Tester = workspace.tester()

#parameters we want to be sampled from
pnames = ['omega_b', 'omega_cdm', 'h', 'tau_reio', 'w0_fld', 'wa_fld', 'N_ur', 'omega_ncdm', 'Omega_k']

#generate cls according to parameter sampling (if not done yet) 
# can take several hours
#Tester.create_cls(pnames,N=120)

#load cls of FULL and NN and create comparisson plots
Tester.compare_cls_full_net(N_lines=1000)
