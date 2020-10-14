import os
import argparse
from ..parameter_sampling import DefaultParamDomain, EllipsoidDomain
from classynet.workspace import GenerationalWorkspace
import classy
import numpy as np

WORKSPACE_DIR = os.path.expanduser("~/CLASSnet_HPC/")

generations = {
    "Net_ST0_Reco":     32,
    "Net_ST0_Reio":     23,
    "Net_ST0_ISW":      23,
    "Net_ST1":          22,
    "Net_ST2_Reco":     24,
    "Net_ST2_Reio":     23,
    "Net_phi_plus_psi": 32,
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("threads", type=int,
                        help="Number of threads."
                             "IMPORTANT: YOU MUST ALSO SET THE ENVIRONMENT VARIABLE OMP_NUM_THREADS"
                             "TO THE SAME VALUE SINCE IT CANNOT BE SET FROM WITHIN PYTHON!"
                        )
    parser.add_argument("-w", "--warmup-iterations", type=int, help="number of warmup iterations", default=2)
    parser.add_argument("-i", "--iterations",        type=int, help="number of actual benchmark iterations", default=25)
    args = parser.parse_args()

    nthreads = args.threads

    workspace = GenerationalWorkspace(WORKSPACE_DIR, generations)

    bm_runner = workspace.benchmark_runner(warmup=args.warmup_iterations, iterations=args.iterations, nthreads=nthreads)
    bm_runner.run()

