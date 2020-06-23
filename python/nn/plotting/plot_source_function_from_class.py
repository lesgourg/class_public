from plot_source_function import plot_source_function
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "../")
import numpy as np

def plot_source_function_from_CLASS(cosmo):
    tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
    sources, k, tau = cosmo.get_sources()
    print(sources.keys())

    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
    plot_source_function(ax[0, 0], k, tau, sources["t0"], tau_rec=tau_rec, title="$S_{T_0}$")
    plot_source_function(ax[0, 1], k, tau, sources["t0_sw"], tau_rec=tau_rec, title="$S_{T_0}^{SW}$")
    plot_source_function(ax[0, 2], k, tau, sources["t0_isw"], tau_rec=tau_rec, title="$S_{T_0}^{ISW}$")
    plot_source_function(ax[1, 0], k, tau, sources["t1"], tau_rec=tau_rec, title="$S_{T_1}$")
    ax[1, 1].remove()
    ax[1, 2].remove()
    plot_source_function(ax[2, 0], k, tau, sources["t2"], tau_rec=tau_rec, title="$S_{T_2}$")
    plot_source_function(ax[2, 1], k, tau, sources["t2_reco"], tau_rec=tau_rec, title="$S_{T_2}^{reco}$")
    plot_source_function(ax[2, 2], k, tau, sources["t2_reio"], tau_rec=tau_rec, title="$S_{T_2}^{reio}$")

    return fig, ax

def plot_source_functions_vs_analytic_approx(cosmo):
    import libbessel
    sources, k, tau = cosmo.get_sources()
    tau_rec = cosmo.get_current_derived_parameters(["tau_rec"])["tau_rec"]
    thermos = cosmo.get_thermos_for_NN()
    tau_th = thermos["tau"]

    r_s, tau_bg = cosmo.get_backgrounds_for_NN()
    r_s = np.interp(tau, tau_bg, r_s)

    cos = np.cos(np.outer(r_s, k))
    sin = np.sin(np.outer(r_s, k))

    j1_n, j2_n = libbessel.compute_bessels(k, tau - tau_rec)

    g = np.interp(tau, np.flip(tau_th), np.flip(thermos["g"]))
    g_reco = np.interp(tau, np.flip(tau_th), np.flip(thermos["g_reco"]))
    g_reio = np.interp(tau, np.flip(tau_th), np.flip(thermos["g_reio"]))

    approx1 = g_reco[:, None] * cos
    approx2 = g_reio[:, None] * j1_n

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))
    plot_source_function(ax[0, 0], k, tau, sources["t2_reco"], tau_rec=tau_rec, title="$S_{T_2}^{reco}$")
    plot_source_function(ax[1, 0], k, tau, approx1.T, tau_rec=tau_rec, title="$g_{reco}(\\tau) \\cos(k r_s)$")
    plot_source_function(ax[0, 1], k, tau, sources["t2_reio"], tau_rec=tau_rec, title="$S_{T_2}^{reio}$")
    plot_source_function(ax[1, 1], k, tau, approx2.T, tau_rec=tau_rec, title="$g_{reio}(\\tau) j_1(k (\\tau - \\tau_{rec}))$")

    return fig, ax

if __name__ == "__main__":
    from classy import Class
    cosmo = Class()

    cosmo_params = {
            "omega_b":   0.02242,
            "omega_cdm": 0.11933,
            "H0":        67.66,
            "tau_reio":  0.0561,
            }

    CLASS_params = {
            'output':'mPk,tCl,pCl,lCl,mTk,vTk,dCl',

        
            "N_ncdm": 1,
            "N_ur": 2.0328,
            "m_ncdm": 0.06,

            "compute damping scale": "yes"
            }

    CLASS_params.update(cosmo_params)
    cosmo.set(CLASS_params)

    cosmo.compute(level=["perturb"])

    fig, _ = plot_source_function_from_CLASS(cosmo)
    # fig.tight_layout()
    # fig.savefig("plots/split_source_function.pdf", bbox_inches="tight")
    fig.savefig("plots/split_source_function.png", dpi=250, bbox_inches="tight")

    fig, _ = plot_source_functions_vs_analytic_approx(cosmo)
    # fig.tight_layout()
    # fig.savefig("plots/analytic_approx.pdf", bbox_inches="tight")
    fig.savefig("plots/analytic_approx.png", dpi=250, bbox_inches="tight")

    plt.show()
