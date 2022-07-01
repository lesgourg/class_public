from classynet.workspace import Workspace
import os

# This script is called by /nn/examples/run_benchmark.sh after running the benchmark.
# All plots will be stored into /path/to/workspace/results/benchmark/

WORKSPACE_DIR = os.path.expanduser("../../../classnet_workspace")

workspace = Workspace(WORKSPACE_DIR)

benchplotter = workspace.benchmark_plotter()

benchplotter.plot_all()
