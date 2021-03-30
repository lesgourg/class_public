import random
import json
import os

import classy

import classynet.utils

class BenchmarkRunner:

    def __init__(self, workspace, nthreads, iterations=25, warmup=5):
        self.workspace = workspace
        self.nthreads = nthreads
        self.iterations = iterations
        self.warmup = warmup

    def run(self):
        settings_warmup = self._load_params(self.warmup)
        settings_actual = self._load_params(self.iterations)

        print("running benchmark for {} threads".format(self.nthreads))

        # warmup
        print("running CLASS {} times (warmup)".format(self.warmup))
        self._run_benchmark(settings_warmup)
        print()
        print("FINISHED WARMUUP")
        print()

        # actual run
        print("running CLASS {} times (benchmark)".format(self.iterations))
        no_nn, nn = self._run_benchmark(settings_actual)
        result = {"nn": nn, "no_nn": no_nn}

        self._save_results(result)

    def _load_params(self, count):
        loader = self.workspace.loader()
        _, cosmo_params = loader.cosmological_parameters()
        cosmo_params = utils.transpose_dict_of_lists(cosmo_params)
        # load fixed parameters to be set for each evaluation
        manifest = loader.manifest()
        fixed = manifest["fixed"]

        result = []
        for p in random.sample(cosmo_params, count):
            p.update(fixed)
            result.append(p)
        return result

    def _run_benchmark(self, list_of_settings):
        if os.environ.get("OMP_NUM_THREADS") != str(self.nthreads):
            raise ValueError("BenchmarkRunner constructed with nthreads={}, but OMP_NUM_THREADS is {} instead!".format(
                self.nthreads, os.environ.get("OMP_NUM_THREADS")))

        def run(use_nn):
            return [self.run_class(settings, use_nn=use_nn) for settings in list_of_settings]

        self._show_msg("performing run WITHOUT neural networks")
        without_nn = run(False)
        print()
        self._show_msg("performing run WITH neural networks")
        with_nn    = run(True)

        return without_nn, with_nn

    def run_class(self, params, use_nn=False, verbose=0):
        cosmo = classy.Class()
        cosmo.set(params)
        if use_nn:
            if verbose>0:
                print("USING NEURAL NETWORKS!")
            cosmo.set({"neural network path": self.workspace})
        timings = {}
        cosmo.compute(level=["perturb"], performance_report=timings)
        cosmo.struct_cleanup()
        return timings

    def _save_results(self, results):
        path = self.workspace.benchmark / "data.json"
        if os.path.exists(path):
            with open(path, "r") as f:
                data = json.load(f)
        else:
            data = {}
        if str(self.nthreads) in data:
            print("WARNING: `{}` already contains data for nthreads={}; overwriting...".format(path, self.nthreads))
        data[str(self.nthreads)] = results
        print("saving benchmark results to", path)
        with open(path, "w") as out:
            json.dump(data, out, indent=4)

    def _show_msg(self, msg):
        print("#" * 80)
        for line in msg.split("\n"):
            print("#", line)
        print("#" * 80)

