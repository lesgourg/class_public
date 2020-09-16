import random
import json
import os

import classy

from . import utils

class BenchmarkRunner:

    def __init__(self, workspace, iterations=25, warmup=5):
        self.workspace = workspace
        self.iterations = iterations
        self.warmup = warmup

    def run(self, thread_counts):
        settings_warmup = self._load_params(self.warmup)
        settings_actual = self._load_params(self.iterations)

        nn = {}
        no_nn = {}

        result = {}
        for thread_count in thread_counts:
            self._show_msg("running benchmark for {} threads".format(thread_count))
            # warmup
            print("running CLASS {} times (warmup)".format(self.warmup))
            self._run_benchmark(thread_count, settings_warmup)
            # actual run
            print("running CLASS {} times (benchmark)".format(self.iterations))
            no_nn, nn = self._run_benchmark(thread_count, settings_actual)
            result[thread_count] = {"nn": nn, "no_nn": no_nn}
        self._save_results(result)

    def _load_params(self, count):
        loader = self.workspace.loader()
        cosmo_params, _ = loader.cosmological_parameters()
        cosmo_params = utils.transpose_dict_of_lists(cosmo_params)
        # load fixed parameters to be set for each evaluation
        manifest = loader.manifest()
        fixed = manifest["fixed"]

        result = []
        for p in random.sample(cosmo_params, count):
            p.update(fixed)
            result.append(p)
        return result

    def _run_benchmark(self, threads, list_of_settings):
        # TODO handle threads
        import importlib

        original = os.environ.get("OMP_NUM_THREADS", None)
        os.environ["OMP_NUM_THREADS"] = "1"

        importlib.reload(classy)

        def run(use_nn):
            return [self.run_class(settings, use_nn=use_nn) for settings in list_of_settings]

        without_nn = run(False)
        with_nn    = run(True)

        if original is not None:
            os.environ["OMP_NUM_THREADS"] = original
        else:
            del os.environ["OMP_NUM_THREADS"]

        return without_nn, with_nn

    def run_class(self, params, use_nn=False):
        cosmo = classy.Class()
        cosmo.set(params)
        if use_nn:
            print("USING NEURAL NETWORKS!")
            cosmo.set({"neural network path": self.workspace})
        timings = {}
        cosmo.compute(level=["perturb"], performance_report=timings)
        cosmo.struct_cleanup()
        return timings

    def _save_results(self, data):
        path = self.workspace.benchmark / "data.json"
        print("saving benchmark results to", path)
        with open(path, "w") as out:
            json.dump(data, out)

    def _show_msg(self, msg):
        print("#" * 80)
        for line in msg.split("\n"):
            print("#", line)
        print("#" * 80)

