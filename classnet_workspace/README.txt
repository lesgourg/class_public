A short introduction to ClassNet,
10.06.2021, Stadtmann


Basics:

This is a ClassNet workspace that can be used to run Class in the Neural Network mode.
In this mode, the Source functions ST0, ST1, ST2, phi+psi, delta_m and delta_cb will be 
estimated with networks if the input parametersa are in a region where the networks are valid.
To do so, pass the parameter "workspace_path" to Class. It should be set to your ClassNet
Workspace directory, like this one: 
"workspace_path":"/path/to/this/class_folder/classnet_workspace/"



Important!!!:

The following files and folders are crutial for ClassNet to work!!! 
More information abotu them below.

models/Net_phi_plus_psi.pt
models/Net_ST0_Reco.pt
models/Net_ST0_Reio.pt
models/Net_ST0_ISW.pt
models/Net_ST1.pt
models/Net_ST2_Reco.pt
models/Net_ST2_Reio.pt
data/k.npy
data/lcdm_11df_sn.bestfit
data/lcdm_11df_sn.covmat
training/normalization.json



Additional Parameters:
(only pass them if you set workspace_path as well)

   nn_verbose
Additionally, you can use "nn_verbose" to toggle the verbose of the ClassNet mode. 
The default is 2 for a big warining pannel and some information about the network evaluation.
For Parameter Estimation or shorter shell output, we recommend using "nn_verbose":1.
If you set nn_verbose to 0, there will be no information whether networks have been used.
If you use this option, make sure that you are not accidentally using the wrong mode.

   nn_cheat
This parameter can be used to calculate components of single source functions with the full Class code.
If you want to calculate for example the reionisation components of ST0 and ST2 with the full code and
all othe rcomponents with the neural networks, you can pass the neural network path along with 
"nn_cheat":set(["t0_reio_no_isw","t2"])
Here is a list of all source function components currently supported by ClassNet:
["t0_reco_no_isw","t0_reio_no_isw","t0_isw","t1","t2_reco","t2_reio","phi_plus_psi","delta_m","delta_cb"]

   network_delta_chi_squared
The neural networks are trained and validated on specific regions of cosmological parameters. 
To determine whether the neural networks can be used for a specific combination of parameters,
a parameter delta_chi_squared is calculated. If this parameter exeeds the validation limit, the
Class code switches to the full calculation. The parameter "network_delta_chi_squared" can be
requested as a derived parameter to obtain information about how far the cosmological parameters 
are inside or outside of the validation region.



The Workspace Structure:

   manifest.json
A file "manifest.json" gives information about the settings passed to Class for the training, 
as well as the cosmological parameters for which the network is trained.
Note that for the training of the current network, we passed additional accuracy parameters for
ClassFull. These were not used for the validation and therefore are not listed int he manifest 
- the tested accuracy referrs to Class run with the default accuracy settings of 3.0.

   models/
The actual trained models are saved in the folder models/. 
When you train the networks, you can specify a generation number to keep track of the
various network versions. However, ClassNet simply uses the networks without this number.
We refer to the neural networks in this workspace at the time of writing this README with
the version number 101, but these will be updated when better networks are found.

   data/k.npy
The folder data/ contains the file k.npy. This file is very important for the usage of ClassNet.
It contains the k array for which the networks have been trained. The file domain.json contains 
information about the domain on which the networks are trained.

   data/lcdm_11df_sn.bestfit and data/lcdm_11df_sn.covmat
The files lcdm_11df_sn.bestfit and lcdm_11df_sn.covmat contain the bestfit and covariance matrix 
of the planck 2018 survey for eleven parameters, nine of them used for the network training
[h, omega_m, omega_cdm, tau_reio, omega_k, w0_fld, wa_fld, N_ur, omega_ncdm]
and two [A_s, n_s] wich are not relevant for the networks. The bestfit and covmat files are used 
to determine the training and validation region and to verify whether a set of cosmological 
parameters is inside this region when executing ClassNet.

   history/
The folder history/ contains the training and validation loss after each epoch of training
and validating. There is one .csv file for each network.

   training/ and validation/
In the training/ folder, there is a file normalization.json. To train the neural networks,
the input and output are normalized. This file contains the information on this normalization,
which is not only neede for training, but also for validation, testing and the usual use of classnet.
Additionally, hese folders contain the parameters for the specific training and validation steps.



