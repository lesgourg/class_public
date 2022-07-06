import os
import argparse
from classynet.workspace import Workspace
import classy
import numpy as np

# This script is called by the bash script /nn/examples/run_benchmark.sh
# It takes the arguments number of warmup rounds and iterations.

WORKSPACE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../classnet_workspace")


parser = argparse.ArgumentParser()
parser.add_argument("threads", type=int,
                    help="Number of threads."
                            "IMPORTANT: YOU MUST ALSO SET THE ENVIRONMENT VARIABLE OMP_NUM_THREADS"
                            "TO THE SAME VALUE SINCE IT CANNOT BE SET FROM WITHIN PYTHON!"
                    )
parser.add_argument("-w", "--warmup-iterations", type=int, help="number of warmup iterations", default=5)
parser.add_argument("-i", "--iterations",        type=int, help="number of actual benchmark iterations", default=50)
args = parser.parse_args()

nthreads = args.threads

workspace = Workspace(WORKSPACE_DIR)

bm_runner = workspace.benchmark_runner(warmup=args.warmup_iterations, iterations=args.iterations, nthreads=nthreads)
bm_runner.run()

# The benchmarks are stored within the workspace directory under /Workspace/results/benchmark