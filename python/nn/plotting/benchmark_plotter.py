import pandas as pd

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np

class BenchmarkPlotter:
    def __init__(self, workspace, plotting_style = 'latex'):
        self.workspace = workspace
        self.print_titles = True
        self.figsize = (6,4)
        if plotting_style == 'latex':
            from matplotlib import rc
            rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 20})
            rc('text', usetex=True)
            self.print_titles = False
        else:
            pass
    
    # plot all 4 performance benchmark plots
    def plot_all(self):
        # load data
        df_no_nn, df_nn = self._load_data()

        self.plot_absolute_perturb_time(df_no_nn, df_nn)
        self.plot_perturb_fraction(df_no_nn, df_nn)
        self.plot_speedup(df_no_nn, df_nn)

        self.figsize = (11,8)
        self.plot_perturb_contributions(df_no_nn, df_nn)

    def plot_absolute_perturb_time(self, df_no_nn, df_nn):
        """
        create bar plots that show the absolute time spent in the perturbation
        module for both Class and ClassNet as a function of the thread number.
        """
        df_perturb = pd.DataFrame({
            "Class":    df_no_nn.perturb,
            "ClassNet": df_nn.perturb
        })

        df_perturb.plot.bar(rot=0)
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.gca().set_yscale("log")
        if self.print_titles:
            plt.title("Time spent in Perturbation Module")
        plt.ylabel("time (s)")
        plt.xlabel("number of threads")
        # NOTE: matplotlib plots vs. the INDEX of the bin, i.e. starting at 0 for threads = 1!
        self._savefig("perturb_abs_time")
        plt.close()

    def plot_perturb_fraction(self, df_no_nn, df_nn):
        """
        create a bar char showing the fraction of the time that is spent inside
        the perturbation module as a function of the number of threads for
        Class and ClassNet.
        """
        df_fraction = pd.DataFrame({
            "Class":    100 * df_no_nn.perturb / df_no_nn.compute,
            "ClassNet": 100 * df_nn.perturb / df_nn.compute,
        })

        df_fraction.plot.bar(rot=0)
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.gca().set_ylim(0,100)
        if self.print_titles:
            plt.title("Fraction of Run Time spent in Perturbation Module")
        plt.ylabel(r"runtime in \texttt{pert} module [\%]")
        plt.xlabel("number of threads")
        plt.legend(loc="upper right")
        self._savefig("perturb_time_fraction")
        plt.close()

    def plot_speedup(self, df_no_nn, df_nn):
        """
        create a bar char showing the speedup of the perturbation module (i.e.
        the ratio `time(Class) / time(ClassNet) - 1`.
        """

        df_speedup = pd.DataFrame({
            "Perturbation": df_no_nn.perturb / df_nn.perturb - 1,
            "Overall EBS": (df_no_nn.compute-df_no_nn.nonlinear) / (df_nn.compute-df_nn.nonlinear) - 1,
        })
        print(df_no_nn.perturb / df_nn.perturb - 1)
        print((df_no_nn.compute-df_no_nn.nonlinear) / (df_nn.compute-df_nn.nonlinear) - 1)
        df_speedup.plot.bar(rot=0)
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.gca().set_yscale("log")
        plt.gca().set_ylabel("speedup")
        plt.gca().set_ylim(1,1000)

        if self.print_titles:
            plt.title("speedup (time(Class) / time(ClassNet) - 1)")
        plt.xlabel("number of threads")
        #plt.gca().get_legend().remove()
        self._savefig("speedup")
        plt.close()

    def plot_perturb_contributions(self, df_no_nn, df_nn):
        def trim_name(name):
            return name[name.find(":")+3:-1]

        # search all individual network contributions
        network_fields = [f for f in df_nn.columns if f.startswith("indiv. network:")]

        # contributions with are displayed in the plot and thus are not required to be considered when calculate the overhead
        subtract_fields = network_fields + [
            "neural network input transformation",
            "neural network output transformation",
            "perturb_init"
        ]
        df_nn["predict overhead"] = df_nn["get all sources"] - df_nn[subtract_fields].sum(axis=1)
        
        # sum CLASSY overhead functions
        classy_overhead = ["update predictor",
        "allocate unused source functions",
        "overwrite source functions",
        "perturb_init",
        "overwrite k array",
        "allocate numpy array of predictions",
            ]
        df_nn["CLASSY overhead"] = df_nn[classy_overhead].sum(axis=1)

        # functions which are connected to the output processing of the NN predictions
        output_comp = ["neural network output transformation",
        "float to double transformation",
            ]
        df_nn["output processing + copying"] = df_nn[output_comp].sum(axis=1)

        fields_start = [
            "neural network input transformation",
        ]
        fields = fields_start + network_fields + ['CLASSY overhead'] + ["output processing + copying"] + ["predict overhead"]
        df_comp = df_nn[fields]

        # sort columns in ascending order of time taken
        sort_key = max(df_comp.index)
        df_comp = df_comp.sort_values(by=sort_key, axis=1, ascending=False)

        # calculate mismatch of summed components to the total computation time due to functions calls or other unconsidered contributions. It is in fact in the magnitude of ~2%
        mismatch = (df_comp.sum(axis=1) - df_nn.perturb) / df_nn.perturb

        # seeked ordering of the contribtions
        my_order = [
        "output processing + copying",  "predict overhead","neural network input transformation", "CLASSY overhead", 
        "indiv. network: 't0_reco_no_isw'",
        "indiv. network: '('phi_plus_psi', 'delta_m', 'delta_cb')'",
        "indiv. network: 't0_isw'",
        "indiv. network: 't1'",
        "indiv. network: 't2_reco'",
        "indiv. network: 't0_reio_no_isw'", "indiv. network: 't2_reio'"]

        df_comp=df_comp.reindex(columns=my_order)
        labels = []

        import re
        pat = re.compile(r"indiv. network: '(.*)'")
        replacements = {
            "neural network input transformation": "input processing",
            "indiv. network: 't0_reco_no_isw'": r"$S_{T_0,\mathrm{reco}}$",
            "indiv. network: 't0_reio_no_isw'": r"$S_{T_0,\mathrm{reio}}$",
            "indiv. network: 't0_isw'": r"$S_{T_0,\mathrm{ISW}}$",
            "indiv. network: 't1'": r"$S_{T_1}$",
            "indiv. network: 't2_reco'": r"$S_{T_2,\mathrm{reco}}$",
            "indiv. network: 't2_reio'": r"$S_{T_2,\mathrm{reio}}$",
            "indiv. network: '('phi_plus_psi', 'delta_m', 'delta_cb')'": r"$S_{\phi+\psi},S_{\delta_m},S_{\delta_{cb}}$",
            "predict overhead": "NN overhead"
        }
        for name in df_comp.columns:
            if name in replacements:
                labels.append(replacements[name])
            else:
                labels.append(name)
        
        # define required colormaps
        cm1 = plt.get_cmap('Reds')
        cm2 = plt.get_cmap('viridis')
        my_colors = [cm1(i) for i in [0.8,0.6,0.4,0.2]] + [cm2(i) for i in [0.2,0.35,0.5,0.65,0.8,0.9,1.0]]

        # create enough unique colors
        df_comp.plot.bar(stacked=True, rot=0, color=my_colors)
        min_key = min(df_comp.index)

        # make some room for legend at top
        plt.gca().set_ylim(0, 2.00 * df_comp.sum(axis=1).loc[min_key])
        plt.gca().legend(labels, ncol=2, mode="expand")
        plt.gca().set_axisbelow(True)
        plt.gca().yaxis.grid(True)
        plt.xlabel("number of threads")
        if self.print_titles:
            plt.title("contributions to Perturbation Module")
        plt.ylabel("time (s)")
        plt.tight_layout()
        self._savefig("time_per_step_and_network")

    def _savefig(self, name, fig=None):
        if fig is None:
            fig = plt.gcf()
            fig.set_size_inches(*self.figsize)

        def save(extension):
            path = self.workspace.benchmark / (name + "." + extension)
            print("saving benchmark plot", path)
            fig.savefig(
                path,
                bbox_inches="tight",
                dpi=100
            )
        save("pdf")

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
            df = df.sort_index()
            return df

        df_no_nn = convert_to_df("no_nn")
        df_nn = convert_to_df("nn")

        return df_no_nn, df_nn





