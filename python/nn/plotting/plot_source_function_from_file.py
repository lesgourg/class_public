from plot_source_function import plot_source_function
import matplotlib.pyplot as plt
import h5py as h5
import sys

if __name__ == "__main__":
    assert len(sys.argv) >= 2
    path = sys.argv[1]

    f = h5.File(path, "r")
    k = f["sampling/k"][()]
    tau = f["sampling/tau"][()]
    tau_rec = f["thermos/tau_rec"][()]

    # convert h5 datasets to numpy arrays for plotting
    sources = {key: dataset[()] for key, dataset in f["sources"].items()}

    fig, ax = plt.subplots(nrows=3, ncols=3)
    plot_source_function(ax[0, 0], k, tau, sources["t0"], tau_rec=tau_rec, title="$S_{T_0}$")
    plot_source_function(ax[0, 1], k, tau, sources["t0_sw"], tau_rec=tau_rec, title="$S_{T_0}^{SW}$")
    plot_source_function(ax[0, 2], k, tau, sources["t0_isw"], tau_rec=tau_rec, title="$S_{T_0}^{ISW}$")
    plot_source_function(ax[1, 0], k, tau, sources["t1"], tau_rec=tau_rec, title="$S_{T_1}$")
    ax[1, 1].remove()
    ax[1, 2].remove()
    plot_source_function(ax[2, 0], k, tau, sources["t2"], tau_rec=tau_rec, title="$S_{T_2}$")
    plot_source_function(ax[2, 1], k, tau, sources["t2_reco"], tau_rec=tau_rec, title="$S_{T_2}^{reco}$")
    plot_source_function(ax[2, 2], k, tau, sources["t2_reio"], tau_rec=tau_rec, title="$S_{T_2}^{reio}$")

    plt.show()

