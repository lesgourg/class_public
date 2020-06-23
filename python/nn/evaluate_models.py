import argparse
import os
import random
from pprint import pprint

import numpy as np
import scipy.interpolate
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import h5py as h5

from tensorflow import keras
import torch

from tabulate import tabulate

import utils
import planck_best_fit
import predictors
import current_transformer
import k_standard
from plotting.plot_source_function import plot_source_function
import common_class_settings
import test_models_on_Cls

import models.T0_reco
from models.original import T0, T1, T2, phi_plus_psi

from models.T0_reco import Net_ST0_Reco
from models.T0_reio import Net_ST0_Reio
from models.T0_isw import Net_ST0_ISW
from models.T1 import Net_ST1
from models.T2_reco import Net_ST2_Reco
from models.T2_reio import Net_ST2_Reio
from models.phi_plus_psi import Net_phi_plus_psi

WORKSPACE_DIR = os.path.expandvars("$CLASSNET_TRAINING_WORKSPACE")

def load_torch_model(path, cls, slicing=None):
    net = cls()
    use_cuda = False
    if use_cuda:
        assert torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    print("moving model to {}".format(device))
    net.load_state_dict(torch.load(path, map_location=device))
    net.to(device)
    net.eval()
    return predictors.TorchModel(net, device, slicing=slicing)

def load_keras_model(path):
    model = keras.models.load_model(path)
    return predictors.KerasModel(model)

class NthAccessor:
    def __init__(self, n):
        self.n = n

    def __call__(self, items):
        return items[0][..., self.n]

class CutT0Reio:
    def __init__(self, threshold=10):
        self.threshold = threshold

    def __call__(self, S, raw_inputs):
        S[:, raw_inputs["tau"] / raw_inputs["tau_rec"] < self.threshold] = 0.
        return S

class CutT0Reco:
    def __init__(self, threshold=2):
        self.threshold = threshold

    def __call__(self, S, raw_inputs):
        S[:, raw_inputs["tau"] / raw_inputs["tau_rec"] > self.threshold] = 0.
        return S

def get_models(choices, new=True):
    ret = get_models_new(choices) if new else get_models_original()

    models = ret[0]

    return ret

def get_models_new(choices):
    wpath = lambda *p: os.path.join(WORKSPACE_DIR, *p)
    ltm = load_torch_model

    reco_slicing = predictors.TimeSlicingReco(4)
    reio_slicing = predictors.TimeSlicingReio(0.6)

    # phi_psi_name = ("phi", "psi", "delta_m", "phi_prime")
    phi_psi_name = ("phi_plus_psi", "delta_m")

    models = {
            "t0_reco_no_isw": ltm(
                wpath("t0_reco",choices["t0_reco_no_isw"], "models", "model.pt"),
                Net_ST0_Reco,
                reco_slicing
            ),
            "t0_reio_no_isw": ltm(
                wpath("t0_reio", choices["t0_reio_no_isw"], "models", "model.pt"),
                Net_ST0_Reio,
                reio_slicing
            ),
            "t0_isw": ltm(
                wpath("t0_isw", choices["t0_isw"], "models", "model.pt"),
                Net_ST0_ISW
            ),
            "t1": ltm(
                wpath("t1", choices["t1"], "models", "model.pt"),
                Net_ST1
            ),
            "t2_reco": ltm(
                wpath("t2_reco", choices["t2_reco"], "models", "model.checkpoint.20.h5"),
                Net_ST2_Reco,
                predictors.TimeSlicingReco(4)
            ),
            "t2_reio": ltm(
                wpath("t2_reio", choices["t2_reio"], "models", "model.pt"),
                Net_ST2_Reio,
                reio_slicing
            ),
            phi_psi_name: ltm(
                wpath("phi_plus_psi", choices["phi_plus_psi"], "models", "model.pt"),
                Net_phi_plus_psi
            ),
            }

    rules = {
            # "t0_isw":           ((phi_psi_name,), compute_t0_isw, True),
            "t0":           (("t0_reco_no_isw", "t0_reio_no_isw", "t0_isw"), sum, False),
            # "t1":           ((("phi", "psi", "delta_m"),), compute_t1, True),
            "t2":           (("t2_reco", "t2_reio"), sum, False),
            "phi_plus_psi": ((phi_psi_name,), NthAccessor(0), False),
            "delta_m":      ((phi_psi_name,), NthAccessor(1), False),
            }

    funcs = {
        # "t0_reco_no_isw": CutT0Reco(threshold=5),
        # "t0_reio_no_isw": CutT0Reio(threshold=10),
        }

    return models, rules, funcs

temp_cosmo_cheat = None

def compute_t1(output, cosmo, tau):
    phi = output[0][..., 0]
    psi = output[0][..., 1]
    psi_minus_phi = psi - phi
    # psi_minus_phi = output[0][..., 1]
    thermo = cosmo.get_thermodynamics()
    tau_th = np.flip(thermo["conf. time [Mpc]"])
    e_kappa_th = np.flip(thermo["exp(-kappa)"])

    e_kappa = np.interp(tau, tau_th, e_kappa_th)

    # t1 = e_kappa[None, :] * k_standard.K_STANDARD[:, None] * (psi - phi)
    t1 = e_kappa[None, :] * k_standard.K_STANDARD[:, None] * psi_minus_phi

    return t1

def compute_t0_isw(output, cosmo, tau):

    # phi = output[0][..., 0]
    # _, phi_prime = np.gradient(phi, k_standard.K_STANDARD, tau)

    t0_isw = output[0][..., 3]
    return t0_isw

    thermo = cosmo.get_thermodynamics()
    tau_th = np.flip(thermo["conf. time [Mpc]"])
    e_kappa_th = np.flip(thermo["exp(-kappa)"])
    e_kappa = np.interp(tau, tau_th, e_kappa_th)

    phi_prime = t0_isw / (2 * e_kappa[None, :])

    # phi_prime = output[0][..., 3]

    src_class, k_class, tau_class = temp_cosmo_cheat.get_sources()
    phi_class = src_class["phi"]
    _, phi_prime_class = np.gradient(phi_class, k_class, tau_class)

    import matplotlib as mpl
    mpl.use("qt4agg")

    def plot_tau_slice():
        k_test = 2e-2
        k_idx_nn = np.argmin(np.abs(k - k_test))
        k_idx_class = np.argmin(np.abs(k_class - k_test))

        plt.figure()
        levels = 200
        plt.subplot(221)
        plot_source_function(plt.gca(), k_standard.K_STANDARD, tau, phi, tau_rec=278, levels=levels)
        plt.subplot(222)
        plot_source_function(plt.gca(), k_standard.K_STANDARD, tau, phi_prime, tau_rec=278, levels=levels)
        plt.subplot(223)
        plot_source_function(plt.gca(), k_class, tau_class, phi_class, tau_rec=278, levels=levels)
        plt.subplot(224)
        plot_source_function(plt.gca(), k_class, tau_class, phi_prime_class, tau_rec=278, levels=levels)
        plt.axvline(k_test, color="k")

        plt.figure(figsize=(10, 10))
        plt.subplot(211)
        plt.title("$\\phi$")
        plt.semilogx(tau, phi[k_idx_nn], label="NN")
        plt.semilogx(tau_class, phi_class[k_idx_class], label="CLASS")
        plt.grid()
        plt.legend()
        plt.subplot(212, sharex=plt.gca())
        plt.title("$\\phi'$")
        plt.semilogx(tau, phi_prime[k_idx_nn], label="NN deriv")
        plt.semilogx(tau_class, phi_prime_class[k_idx_class], label="CLASS deriv")
        plt.xlabel(r"$\tau$")
        plt.grid()
        plt.legend()

        plt.twinx()
        plt.hist(tau, bins=np.geomspace(tau[0], tau[-1], 75), color="g", alpha=0.5, label="$\\tau$ sampling histogram")
        plt.legend()

        plt.savefig("../plots/phi_prime_from_phi.png", dpi=250)
        plt.show()

    def plot_k_slice():
        tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
        tau_rec_idx = np.argmin(np.abs(tau_rec - tau))
        tau_idx = tau_rec_idx + 30

        fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(6, 8))
        axes[0].set_title("T0_isw network prediction")
        plot_source_function(axes[0], k_standard.K_STANDARD, tau, t0_isw, tau_rec=tau_rec, show_colorbar=False)
        axes[0].axhline(tau[tau_idx] / tau_rec, color="k", alpha=0.5)
        axes[0].axvline(2e-2, color="k", ls=":", alpha=0.5)

        axes[1].set_title(r"$k \phi'_{{CLASS}}$ @ $\tau / \tau_{{rec}} = {:.3f}$".format(tau[tau_idx] / tau_rec))
        axes[1].semilogx(k_class, k_class * phi_prime_class[:, tau_idx])
        axes[1].axvline(2e-2, color="k", ls=":", alpha=0.5)
        axes[1].grid()
        axes[1].set_xlabel("$k$")

        hist_ax = axes[1].twinx()
        hist_ax.hist(k_standard.K_STANDARD,
                     bins=np.geomspace(k_standard.K_STANDARD[0], k_standard.K_STANDARD[-1], 50),
                     color="g", alpha=0.5, label="$k$ sampling histogram")
        hist_ax.legend()

        plt.show()

    # plot_tau_slice()
    # plot_k_slice()

    return t0_isw

def extract_phi_plus_psi(output):
    phi = output[0][..., 0]
    psi = output[0][..., 1]
    return phi + psi
    # psi_minus_phi = output[0][..., 1]
    # phi+psi = 2 * phi + (psi - phi)
    # return 2 * phi + psi_minus_phi

def extract_delta_m(output):
    return output[0][..., 2]

def get_models_original():
    wpath = lambda *p: os.path.join(WORKSPACE_DIR, *p)
    ltm = load_torch_model

    models = {
            "t0":                    predictors.TorchModel(T0.Net_ST0(), torch.device("cpu")),
            "t1":                    predictors.TorchModel(T1.Net_ST1(), torch.device("cpu")),
            "t2":                    predictors.TorchModel(T2.Net_ST2(), torch.device("cpu")),
            ("phi_plus_psi", "delta_m"): predictors.TorchModel(phi_plus_psi.Net_phi_plus_psi(), torch.device("cpu")),
            }

    rules = {
            "phi_plus_psi": ((("phi_plus_psi", "delta_m"),), NthAccessor(0)),
            "delta_m":      ((("phi_plus_psi", "delta_m"),), NthAccessor(1)),
            }

    funcs = {
        }

    return models, rules, funcs


def get_predictor_factory(choices, cheat=False, cosmo_cheat=None, cheat_on=None):
    input_transformer = current_transformer.get_input_transformer_normalizer()
    target_transformer = current_transformer.get_target_transformer_normalizer(k=k_standard.K_STANDARD)

    models, rules, funcs = get_models(choices)

    def predictor_factory(cosmo_params):
        predictor = predictors.TreePredictor(
                cosmo_params,
                input_transformer=input_transformer,
                target_transformer=target_transformer,
                models=models,
                rules=rules,
                funcs=funcs,
                )

        if cheat:
            predictor.cheat = cheat
            predictor.cosmo_cheat = cosmo_cheat
            predictor.cheat_on = cheat_on

        return predictor

    return predictor_factory

def get_testing_parameters(i):
    """
    Return the cosmological parameters of the `i`th cosmology in the testing set.
    """
    with h5.File(os.path.expandvars("$CLASSNET_DATA/cosmological_parameters_lhs_testing.h5"), "r") as f:
        cosmo_params = {name: float(arr[i]) for name, arr in f.items()}
    return cosmo_params

def get_cosmological_parameters():
    """
    Return the cosmological parameters of the cosmology for which to produce
    the plots.
    """
    return get_testing_parameters(42)
    # return planck_best_fit.planck_best_fit_mean

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str, help="output directory for plots, etc.")
    parser.add_argument("--models", type=str, help="list of models and model ids to evaluate; format: `model_name_1:model_id_1;model_name_2:model_id_2;<...>")
    parser.add_argument("-s", "--summary", action="store_true",
            help="Produce a plot containing all source functions (both from CLASS and networks + residual)")
    parser.add_argument("--test", action="store_true",
            help="will perform test of Cls")
    parser.add_argument("--test-on", type=int, default=100,
            help="number of cosmologies from testing set to test on")
    parser.add_argument("--cheat", nargs="+", help="quantities for which not to use NN, but 'cheat' using CLASS")
    parser.add_argument("--skip-summary", action="store_true", help="Skip the plots of the source functions")
    return parser

def run_class(cosmo_params):
    cosmo = common_class_settings.get_instance()
    cosmo.set(cosmo_params)
    cosmo.set(lensing="yes")
    cosmo.compute()
    return cosmo

def run_class_nn(cosmo_params, predictor_factory):
    cosmo = common_class_settings.get_instance()
    cosmo.set(cosmo_params)
    # cosmo.set(lensing="yes")
    predictor = predictor_factory(cosmo_params)
    cosmo.enable_NN(predictor)
    perf = {}
    cosmo.compute(performance_report=perf)

    rows = sorted(perf.items(), key=lambda item: item[1], reverse=True)
    print(tabulate(rows, headers=["step", "time (s)"], tablefmt="fancy_grid"))

    return cosmo

def compute_residual(k, tau, S, k_, tau_, S_):
    assert S.shape == (k.size, tau.size)
    spline = scipy.interpolate.RectBivariateSpline(k_, tau_, S_)
    return spline(k, tau) - S

def create_summary_plots(k, tau, sources, k_nn, tau_nn, sources_nn, tau_rec=None, labels=None, k_cut=0.42):
    n_plots = len(sources)

    if labels is None:
        labels = [None] * n_plots

    fig, axes = plt.subplots(nrows=3, ncols=n_plots, figsize=(1.2 * n_plots * 4, 3 * 4), squeeze=False)

    k_mask = k < k_cut
    k_nn_mask = k_nn < k_cut
    k = k[k_mask]
    k_nn = k_nn[k_nn_mask]

    for name, label, axes_col in zip(sources, labels, axes.T):
        print(f"Plotting source function '{name}'")
        S = sources[name]
        S_ = sources_nn[name]

        S = S[k_mask]
        S_ = S_[k_nn_mask]

        residual = compute_residual(k, tau, S, k_nn, tau_nn, S_)
        rms = np.sqrt(np.mean(np.square(residual)))

        levels = 256
        plot_source_function(axes_col[0], k, tau, S, tau_rec=tau_rec, levels=levels, title=label + " CLASS")
        plot_source_function(axes_col[1], k_nn, tau_nn, S_, tau_rec=tau_rec, levels=levels, title=label + " NN")
        plot_source_function(axes_col[2], k, tau, residual, tau_rec=tau_rec, levels=levels, title=label + " diff")
        # axes_col[2].text(.95, .05, "RMS: {:.2e}".format(rms), transform=axes_col[2].transAxes, va="bottom", ha="right", alpha=0.5)

    for ax in axes.flat:
        ax.label_outer()

    return fig

# TODO this is essentially a duplicate of the above
def create_source_function_plot(k, tau, source, k_nn, tau_nn, source_nn, tau_rec=None, label=None, k_cut=0.42):
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(3 * 4, 4))

    k_mask = k < k_cut
    k_nn_mask = k_nn < k_cut
    k = k[k_mask]
    k_nn = k_nn[k_nn_mask]

    print(f"Plotting source function '{label}'")
    S = source
    S_ = source_nn

    S = S[k_mask]
    S_ = S_[k_nn_mask]

    residual = compute_residual(k, tau, S, k_nn, tau_nn, S_)
    rms = np.sqrt(np.mean(np.square(residual)))

    levels = 256
    plot_source_function(axes[0], k, tau, S, tau_rec=tau_rec, levels=levels, title=label + " CLASS")
    plot_source_function(axes[1], k_nn, tau_nn, S_, tau_rec=tau_rec, levels=levels, title=label + " NN")
    plot_source_function(axes[2], k, tau, residual, tau_rec=tau_rec, levels=levels, title=label + " diff")
    axes[2].text(.95, .05, "RMS: {:.2e}".format(rms), transform=axes[2].transAxes, va="bottom", ha="right", alpha=0.7)

    for ax in axes.flat:
        ax.label_outer()

    return fig


def cut_low_l(cls):
    return {key: value[2:] for key, value in cls.items()}

def create_cl_plot(cls, cls_nn, fields=("tt",)):
    figsize = (len(fields) * 6, 2 * 3)
    fig, axes = plt.subplots(nrows=2, ncols=len(fields), squeeze=False, figsize=figsize)
    print("hi")
    fig.tight_layout()

    cls = cut_low_l(cls)
    cls_nn = cut_low_l(cls_nn)

    ell = cls["ell"]

    factor = ell * (ell + 1) / (2 * np.pi)
    xlabel = "$\\ell$"
    ylabel_prefix = "$\\ell (\\ell + 1) / \\sqrt{2 \\pi}$"

    for field, ax_col in zip(fields, axes.T):
        cl_label = "$C_\ell^{{{}}}$".format(field.upper())
        ylabel = ylabel_prefix + cl_label

        ax_col[0].set_title(cl_label)
        ax_col[0].set_ylabel(ylabel)

        ax_col[0].semilogx(ell, factor * cls[field], color="b", alpha=0.7, label="CLASS")
        ax_col[0].semilogx(ell, factor * cls_nn[field], color="r", alpha=0.7, label="NN")
        ax_col[0].legend()

        if field == "te":
            residual = factor * (cls_nn[field] - cls[field])
            error_label = "$\\ell (\ell + 1) / \\sqrt{2\pi} \\Delta C_\\ell$"
        else:
            residual = (cls_nn[field] - cls[field]) / cls[field]
            error_label = "$\\Delta C_\\ell / C_\\ell^{ref}$"

        ax_col[1].semilogx(ell, residual, color="r", alpha=0.7, label="residual")
        ax_col[1].axhline(0, color="b", alpha=0.7)
        ax_col[1].legend()
        ax_col[1].grid()

        ax_col[1].set_xlabel(xlabel)
        ax_col[1].set_ylabel(error_label)

        xscale = "log" if field == "tt" else "linear"
        for cell in ax_col.flat:
            cell.set_xscale(xscale)

    return fig

def create_cl_plot_from(cosmo, cosmo_nn, fields=("tt",), lensed=False, l_max=2500):
    if not lensed:
        cls = cosmo.raw_cl(l_max)
        cls_nn = cosmo_nn.raw_cl(l_max)
    else:
        cls = cosmo.lensed_cl(l_max)
        cls_nn = cosmo_nn.lensed_cl(l_max)
    return create_cl_plot(cls, cls_nn, fields=fields)

def get_sources_from_predictor(predictor_factory, names, cosmo_nn, cosmo_params, tau, add_k0=False):
    predictor = predictor_factory(cosmo_params)
    k_nn = predictor.get_k(add_k0=add_k0)
    return {name: predictor.predict(name, cosmo_nn, tau, add_k0=add_k0) for name in names}, k_nn, tau

def get_sources_from_class(cosmo, names):
    sources, _, _ = cosmo.get_sources()
    return {name: sources[name] for name in names}

def perform_Cl_test(count, predictor_factory, output_dir):
    test_indices = random.sample(range(1000), count)
    cosmo_params = list(map(get_testing_parameters, test_indices))
    print("Testing networks on {} cosmologies from testing set.".format(len(cosmo_params)))
    (cosmo_params, Cls, Cls_nn, errs), fig_cl_error, figs_cl_error, fig_pk = test_models_on_Cls.test_models(
            cosmo_params,
            test_indices,
            predictor_factory,
            select=["tt", "te", "ee"],
            use_cache=False,
            )
    # fig_cls.savefig(os.path.join(plot_dir, "Cls_testing.pdf"), bbox_inches="tight")

    with h5.File(os.path.join(output_dir, "Cl_error.h5"), "w") as out:
        # Write cosmo params
        g_params = out.create_group("cosmological_parameters")
        param_names = cosmo_params[0].keys()
        for name in param_names:
            arr = np.array([params[name] for params in cosmo_params])
            g_params.create_dataset(name, data=arr)

        errors = errs["tt"]

        out.create_dataset("errors", data=errors)

    return fig_cl_error, figs_cl_error, fig_pk

def plot_mPk(k, cosmo, cosmo_nn):
    pk = np.vectorize(cosmo.pk)(k, 0)
    vec_nn = np.vectorize(lambda k: cosmo_nn.pk(k, 0.0))
    pk_nn = vec_nn(k)

    fig, axes = plt.subplots(nrows=2, sharex=True)

    axes[0].semilogx(k, pk, "b-", label="CLASS", alpha=0.7)
    axes[0].semilogx(k, pk_nn, "r-", label="NN", alpha=0.7)
    axes[0].set_ylabel("$P(k)$")
    axes[0].set_yscale("log")
    axes[0].legend()

    axes[1].axhline(0, color="b")
    axes[1].semilogx(k, (pk_nn - pk) / pk, "r-", label="NN", alpha=0.7)
    axes[1].set_ylabel("$(P_{NN}(k) - P_{CLASS}(k)) / P_{CLASS}(k)$")
    # axes[1].set_yscale("symlog", lintreshy=1e-2)
    axes[1].set_xlabel("$k$ [Mpc${}^{-1}$]")

    fig.tight_layout()

    return fig

if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    print(args)

    # Create output directory if it does not exist yet
    if not os.path.exists(args.output_dir):
        print(f"Creating {args.output_dir} since it does not exist yet.")
        os.makedirs(args.output_dir)

    cosmo_params = get_cosmological_parameters()

    ############################################################
    # 1) compute & plot Cl's
    ############################################################

    with utils.timing("Running CLASS (full)"):
        cosmo = run_class(cosmo_params)

    temp_cosmo_cheat = cosmo

    # This is not nice, but quick and dirty
    sources_class, k_class, tau_class = cosmo.get_sources()

    models.T0_reco.RESULT = sources_class["t0_reco_no_isw"].T
    models.T0_reco.RESULT_K = k_class

    models.T0_isw.RESULT = sources_class["t0_isw"].T
    models.T0_isw.RESULT_K = k_class

    models.phi_plus_psi.RESULT = np.stack((sources_class["phi_plus_psi"].T, sources_class["delta_m"].T), axis=2)
    models.phi_plus_psi.RESULT_K = k_class

    choices = dict(map(tuple, (entry.split(":") for entry in args.models.split(";"))))
    print("CHOSEN MODELS:", choices)

    if args.cheat:
        cheat_on = args.cheat
        print("CHEATING ON", cheat_on)
        predictor_factory = get_predictor_factory(choices, cheat=True, cosmo_cheat=cosmo, cheat_on=cheat_on)
    else:
        predictor_factory = get_predictor_factory(choices)

    sources, k, tau = cosmo.get_sources()

    with utils.timing("Running CLASS (NN)"):
        # Warmup run
        # cosmo_nn = run_class_nn(cosmo_params, predictor_factory)
        cosmo_nn = run_class_nn(cosmo_params, predictor_factory)

    fig_cl = create_cl_plot_from(cosmo, cosmo_nn, fields=("tt", "ee", "te"), lensed=False)
    fig_cl.savefig(os.path.join(args.output_dir, "Cls.pdf"), bbox_inches="tight")

    ############################################################
    # 1.5) if requested, perform test on Cls from testing set
    ############################################################
    if args.test:
        with utils.timing("Testing on Cls..."):
            fig_cl_error, figs_cl_error, fig_pk = perform_Cl_test(args.test_on, predictor_factory, args.output_dir)
            fig_cl_error.savefig(os.path.join(args.output_dir, "Cls_testing_error.pdf"), bbox_inches="tight")
            fig_pk.savefig(os.path.join(args.output_dir, "mPk_error.pdf"), bbox_inches="tight")
            fig_pk.savefig(os.path.join(args.output_dir, "mPk_error.png"), bbox_inches="tight", dpi=200)
            for name, fig in zip(["TT", "EE", "TE"], figs_cl_error):
                fig.savefig(os.path.join(args.output_dir, f"Cls_{name}.pdf"), bbox_inches="tight")
                fig.savefig(os.path.join(args.output_dir, f"Cls_{name}.png"), bbox_inches="tight", dpi=200)

    ############################################################
    # 2) plot source functions
    ############################################################
    if not args.skip_summary:
        names = ["t0_reco_no_isw", "t0_reio_no_isw", "t0_isw", "t1", "t2_reco", "t2_reio", "phi_plus_psi", "delta_m"]

        sources, k, tau = cosmo.get_sources()
        sources = {k: sources[k] for k in names}
        sources_nn, k_nn, tau_nn = get_sources_from_predictor(predictor_factory, names, cosmo_nn, cosmo_params, tau)

        tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]

        with utils.timing("Creating summary plots..."):
            fig_summary = create_summary_plots(k, tau, sources, k_nn, tau_nn, sources_nn, tau_rec=tau_rec, labels=names)
        fig_summary.savefig(os.path.join(args.output_dir, "sources.png"), dpi=250, bbox_inches="tight")

        # Plot the source functions individually
        with utils.timing("Creating individual source function plots..."):
            for field in names:
                fig = create_source_function_plot(
                    k, tau, sources[field],
                    k_nn, tau_nn, sources_nn[field],
                    tau_rec=tau_rec, label=field,
                )
                fig.savefig(os.path.join(args.output_dir, f"source_{field}.png"), dpi=250, bbox_inches="tight")


    ############################################################
    # 2) plot mPk
    ############################################################
    with utils.timing("Plotting P(k)..."):
        # k_pk = np.geomspace(k_standard.K_STANDARD[0], k_standard.K_STANDARD[-1], 500, endpoint=True)
        k_pk = k_standard.K_STANDARD
        plot_mPk(k_pk, cosmo, cosmo_nn).savefig(os.path.join(args.output_dir, "mPk.pdf"))

