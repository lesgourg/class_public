from cobaya.run import run
from cobaya.log import LoggedError
from mpi4py import MPI
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#
# This script showcases how to use CLASSnet within the cobaya framework to generate a MCMC.
# The approach is straight forward by providing the extra parameters such as the workspace path and Setting the use_nn flag to 'yes'.
# The test example explores the parameter space of the LambdaCDM model using the Planck temperature and polarization likelihoods.
# 
# For parallelized way run with: mpiexec -n 8 -x OMP_NUM_THREADS=2 python example_mcmc
#
cobaya_info = {
    "debug": False,
    "output": "./my_mcmc/",
    "timing":True,
    "theory": {
        "classy": {
            "extra_args": {
                # default CLASS settings:
                "output": "tCl,pCl,lCl,mPk",
                "non linear": "halofit",
                "lensing":"yes",
                "P_k_max_1/Mpc":100,
                "l_max_scalars":3000,
                "N_ncdm": 1,
                "compute damping scale":"yes",
                "deg_ncdm":3,
                "N_ur": 0.00641,
                "omega_ncdm": 0.06/93.1,

                # additional CLASSnet settings:
                "use_nn": "yes",
                "workspace_path": os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../classnet_workspace"), 
                "nn_verbose": 3               
                },
            "path": "../../.."
            },
        },
    "likelihood": {
        "planck_2018_lowl.TT":{},
        "planck_2018_lowl.EE":{},
        "planck_2018_highl_plik.TTTEEE":{},
        "planck_2018_lensing.clik":{},
        },
    "params": {
        "Omega_Lambda":{
            "value": 0.,
            "latex": r"\Omega_\Lambda",
            },
        "A_s": {
            "value":"lambda logA: 1e-10*np.exp(logA)",
            "latex":r"A_\mathrm{s}",
            },
        "n_s": {
            "prior": {
                "min": 0.8,
                "max": 1.2,
                },
            "ref": {
                "dist":"norm",
                "loc":0.965,
                "scale":0.004,
                },
            "proposal":0.002,
            "latex":r"n_\mathrm{s}"
            },
        # "H0":{
        #     "value:":"lambda h:h*100.0",
        #     "latex": r"H_0",
        #     },
        "h":{
            "prior":{
                "min":0.20,
                "max":1.00,
                },
            "ref":{
                "dist":"norm",
                "loc":0.67,
                "scale":0.02,
                },
            "proposal":0.01,
            "latex":r"h",
            },
        
        "omega_b":{
            "prior": {
                "min":0.005,
                "max":0.1,
                },
            "ref":{
                "dist":"norm",
                "loc":0.0224,
                "scale":0.0001,
                },
            "proposal":0.0001,
            "latex":r"\Omega_\mathrm{b} h^2"
            },
        "omega_cdm": {
            "prior": {
                "min":0.001,
                "max":0.99,
                },
            "ref": {
                "dist":"norm",
                "loc":0.12,
                "scale":0.001,
                },
            "proposal": 0.0005,
            "latex":r"\Omega_\mathrm{c} h^2",
            },
        "tau_reio": {
            "prior":{
                "min": 0.01,
                "max": 0.8,
                },
            "ref":{
                "dist":"norm",
                "loc":0.055,
                "scale":0.006,
                },
            "proposal":0.003,
            "latex":r"\tau_\mathrm{reio}"
            },
        "logA": {
            "prior": {
                "min":1.61,
                "max":3.91,
                },
            "ref": {
                "dist": "norm",
                "loc":3.05,
                "scale":0.001,
                },
            "proposal":0.001,
            "latex":r"\log(10^{10} A_\mathrm{s})",
            "drop":True
            },
        
        "Omega_m":{
                "latex":r"\Omega_\mathrm{m}",
                },
        "nn_chi2":{
                "latex":r"\triangle_\chi^2",
                },
        },
    "sampler": { 
        "mcmc": {
            "covmat":"auto",
            "drag":True,
            "oversample_power": 0.4,
            "proposal_scale":1.9,
            "Rminus1_stop": 0.05,
            },
        },
    }

success=False
try:
    info=cobaya_info
    updated_info, sampler = run(info)
    success = True
except LoggedError as err:
    pass
success = all(comm.allgather(success))
if not success and rank ==0:
    print("Sampling failed!")
