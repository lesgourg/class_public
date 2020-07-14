import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

def load_errors(path):
    with open(path, "rb") as f:
        return pickle.load(f)


PRETTY_LABELS = {
    "N_ur": r"$N_{ur}$",
    "Omega_k": r"$\Omega_k$",
    "h": r"$h$",
    "omega_b": r"$\omega_b$",
    "omega_cdm": r"$\omega_\mathrm{cdm}$",
    "omega_ncdm": r"$\omega_\mathrm{ncdm}$",
    "tau_reio": r"$\tau_\mathrm{reio}$",
    "w0_fld": r"$w_0$",
    "wa_fld": r"$w_a$",
}

def get_pretty_label(name):
    return PRETTY_LABELS.get(name, name)

def triangle_plot(param_names, param_values, errors, rasterized=True):
    """
    param_names: List[str]
    param_values: List[List[float]] <- inner list must have same order as param_names
    errors: List[float]
    """

    size = len(param_names)
    fig, axes = plt.subplots(nrows=size, ncols=size, figsize=(3 * size, 3 * size))
    fig.subplots_adjust(wspace=0, hspace=0)

    err_mean = np.mean(errors)
    err_std = np.std(errors)

    # err_color = np.clip(errors, err_mean - 3 * err_std, err_mean + 3 * err_std)
    err_color = errors

    # remove plots above diagonal
    for i, j in zip(*np.triu_indices_from(axes, 1)):
        axes[i, j].set_visible(False)

    for i, j in zip(*np.tril_indices_from(axes)):
        ax = axes[i, j]
        if i == size - 1:
            param_x = get_pretty_label(param_names[j])
            ax.set_xlabel(param_x)
        if j == 0:
            param_y = get_pretty_label(param_names[i])
            ax.set_ylabel(param_y)
        ax.label_outer()

        if i == j:
            x = param_values[:, j]
            ax_t = ax.twinx()
            ax_t.scatter(x, errors, c=err_color, rasterized=rasterized)
            ax.set_xlim(x.min(), x.max())
            if i == 0:
                ax.set_ylim(x.min(), x.max())
            else:
                ax.get_yaxis().set_ticks([])
        else:
            x = param_values[:, j]
            y = param_values[:, i]
            ax.scatter(x, y, c=err_color, rasterized=rasterized)
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(y.min(), y.max())

    return fig


def get_plotter(param_names, param_values):
    def plot(qty, title, path):
        print("creating plot for ", title)
        fig = triangle_plot(param_names, param_values, qty)
        fig.suptitle(title)
        print("writing plot to", path)
        fig.savefig(path, bbox_inches="tight")
    return plot

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str, help=".pickle file containing cl and pk errors")
    args = parser.parse_args()
    input_path = Path(args.filename)
    assert input_path.exists()
    assert input_path.is_file()

    output_dir = input_path.parent


    """
        stat_dict = {
            "parameters": cosmo_params,
            "cl_error": cl_err,
            "pk_error": pk_err,
            "cl_error_relative": cl_err_relative,
            "pk_error_relative": pk_err_relative
        }
    """

    ell = np.arange(2, 3000 + 1)
    cv = np.sqrt(2 / (2 * ell + 1))

    errors = load_errors(args.filename)
    params0 = errors[0]["parameters"]
    param_names = list(sorted(params0))

    param_values = np.array([[d["parameters"][key] for key in param_names] for d in errors])

    plot = get_plotter(param_names, param_values)
    plot(
        np.log10(np.array([np.sqrt(np.mean(np.square(d["cl_error_relative"]["tt"]))) for d in errors])),
        title="RMS relative error of $C_\ell^{TT}$ (log10)",
        path=output_dir / "cl_tt_rel_err_rms.pdf"
    )

    plot(
        np.log10(np.array([np.sqrt(np.mean(np.square(d["cl_error_relative"]["tt"] / cv))) for d in errors])),
        title="RMS relative error relative to cosmic variance of $C_\ell^{TT}$ (log10)",
        path=output_dir / "cl_tt_rel_err_cv_rms.pdf"
    )

    plot(
        np.log10(np.array([np.max(np.abs(d["cl_error_relative"]["tt"])) for d in errors])),
        title="max. relative error (over all $\\ell$) of $C_\\ell^{TT}$ (log10)",
        path=output_dir / "cl_tt_rel_err_max.pdf"
    )

    plot(
        np.log10(np.array([np.max(np.abs(d["cl_error_relative"]["tt"] / cv)) for d in errors])),
        title="max. relative error relative to cosmic variance (over all $\\ell$) of $C_\\ell^{TT}$ (log10)",
        path=output_dir / "cl_tt_rel_err_cv_max.pdf"
    )

    plot(
        np.log10(np.array([np.sqrt(np.mean(np.square(d["pk_error_relative"]))) for d in errors])),
        title="RMS of relative error (over all $k$) of $P(k)$ (log10)",
        path=output_dir / "pk_rel_err_rms.pdf"
    )

    plot(
        np.log10(np.array([np.max(np.abs(d["pk_error_relative"])) for d in errors])),
        title="max. relative error (over all $k$) of $P(k)$ (log10)",
        path=output_dir / "pk_rel_err_max.pdf"
    )

