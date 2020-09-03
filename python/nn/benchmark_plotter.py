import pandas as pd

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

class BenchmarkPlotter:
    def __init__(self, workspace):
        self.workspace = workspace

    def plot_all(self):
        pass

    def plot_perturbation_module(self):
        df_no_nn, df_nn = self._load_data()

        plot_functions = [
            self.plot_absolute_perturb_time,
            self.plot_perturb_fraction,
            self.plot_perturb_speedup,
            self.plot_perturb_contributions,
        ]

        for pf in plot_functions:
            pf(df_no_nn, df_nn)

    def plot_absolute_perturb_time(self, df_no_nn, df_nn):
        """
        create bar plots that show the absolute time spent in the perturbation
        module for both CLASS and CLASSnet as a function of the thread number.
        """
        df_perturb = pd.DataFrame({
            "CLASS":    df_no_nn.perturb,
            "CLASSnet": df_nn.perturb
        })

        df_perturb.plot.bar(rot=0)
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.title("Time spent in Perturbation Module")
        plt.ylabel("time (s)")
        plt.xlabel("number of threads")
        self._savefig("perturb_abs_time")
        plt.close()

    def plot_perturb_fraction(self, df_no_nn, df_nn):
        """
        create a bar char showing the fraction of the time that is spent inside
        the perturbation module as a function of the number of threads for
        CLASS and CLASSnet.
        """
        df_fraction = pd.DataFrame({
            "CLASS":    100 * df_no_nn.perturb / df_no_nn.compute,
            "CLASSnet": 100 * df_nn.perturb / df_nn.compute,
        })

        df_fraction.plot.bar(rot=0)
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.title("Fraction of Run Time spent in Perturbation Module")
        plt.ylabel("%")
        plt.xlabel("number of threads")
        self._savefig("perturb_time_fraction")
        plt.close()

    def plot_perturb_speedup(self, df_no_nn, df_nn):
        """
        create a bar char showing the speedup of the perturbation module (i.e.
        the ratio `time(CLASS) / time(CLASSnet) - 1`.
        """
        df_speedup = pd.DataFrame({
            "CLASSnet": df_no_nn.perturb / df_nn.perturb - 1,
        })

        df_speedup.plot.bar(rot=0)
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.title("Speedup (time(CLASS) / time(CLASSnet) - 1)")
        plt.ylabel("time (s)")
        plt.xlabel("number of threads")
        plt.gca().get_legend().remove()
        self._savefig("perturb_speedup")
        plt.close()

    def plot_perturb_contributions(self, df_no_nn, df_nn):
        def trim_name(name):
            return name[name.find(":")+3:-1]

        network_fields = [f for f in df_nn.columns if f.startswith("indiv. network:")]
        fields_start = [
            "perturb_init",
            "neural network input transformation",
        ]
        fields_end = ["neural network output transformation"]
        fields = fields_start + network_fields + fields_end

        labels = fields_start + [trim_name(name) for name in network_fields] + fields_end

        df_comp = df_nn[fields]
        # sort columns in ascending order of time taken (for row corresponding to 1 thread)
        df_comp = df_comp.sort_values(by=1, axis=1, ascending=False)

        labels = []
        import re
        pat = re.compile(r"indiv. network: '(.*)'")
        replacements = {
            "neural network output transformation": "output processing",
            "neural network input transformation": "input processing",
        }
        for name in df_comp.columns:
            match = pat.match(name)
            if match:
                labels.append(match.group(1))
            elif name in replacements:
                labels.append(replacements[name])
            else:
                labels.append(name)

        df_comp.plot.bar(stacked=True, rot=0)
        plt.gca().set_ylim(0, 1.5 * df_comp.sum(axis=1).loc[1])
        plt.gca().legend(labels, ncol=2, mode="expand")
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.xlabel("Number of threads (OMP_NUM_THREADS)")
        plt.title("Contributions to Perturbation Module")
        plt.ylabel("Time (s)")
        plt.tight_layout()
        self._savefig("time_per_step_and_network")

    def _savefig(self, name, fig=None):
        if fig is None:
            fig = plt.gcf()
        path = self.workspace.benchmark / (name + ".png")
        print("saving benchmark plot", path)
        fig.savefig(
            path,
            bbox_inches="tight",
            # TODO increase for final
            dpi=100
        )

    def _load_data(self):
        # dict with format `data[thread_count]["nn" | "no_nn"]
        data = self.workspace.loader().benchmark_data()

        # create a dataframe from the loaded data which has the following shape:
        #          background   compute     input   perturb  perturb_init  thermodynamics
        # threads
        # 1          0.049667  2.609824  0.000248  2.429060      2.429051        0.130781
        # 4          0.049778  2.597886  0.000260  2.417113      2.417105        0.130652
        # ...

        def convert_to_df(name):
            df = pd.DataFrame()
            for thread_count, samples in data.items():
                mean = pd.DataFrame(samples[name]).mean()
                mean.name = int(thread_count)
                df = df.append(mean)
            df.index.name = "threads"
            return df

        df_no_nn = convert_to_df("no_nn")
        df_nn = convert_to_df("nn")
        print(df_nn)

        return df_no_nn, df_nn





