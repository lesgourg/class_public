import numpy as np
import matplotlib.pyplot as plt

def plot_Cls(ref, nn, field="tt", errors=None, residuals=False):
    """
    Plot Cls given a reference (i.e. full computation) `ref` (a dictionary
    obtained from a CLASS instance using e.g. `.raw_cl()`) and a list/generator
    yielding Cls obtained using neural networks.
    `fields` can either be a list of which Cls to plot or just a single value.
    """
    if nn is None:
        nn = []
    elif not isinstance(nn, list):
        nn = [nn]

    fig, axes = plt.subplots(1 + int(residuals), sharex=True)

    ax = axes[0]

    if not residuals:
        ax.set_xlabel(r"$\ell$")

    ax.set_ylabel(r"$\ell (\ell + 1) C_\ell^{{{}}} / (2 \pi)$".format(field.upper()))
    if errors:
        ax.fill_between(ref["ell"][2:], ref[field][2:] - errors[2:], ref[field][2:] + errors[2:], facecolor="k", alpha=0.2)

    l_ref = ref["ell"]
    to_y = lambda c: c * (l_ref * (l_ref + 1)) / (2 * np.pi)

    y_ref = to_y(ref[field])
    ax.semilogx(l_ref[2:], y_ref[2:], label="CLASS")
    for Cls in nn:
        l = Cls["ell"]
        ax.semilogx(l[2:], to_y(Cls[field])[2:], ls="--", color="r", label="CLASSnet")
    ax.legend()

    if residuals:
        ax = axes[1]
        ax.set_xlabel(r"$\ell$")
        ax.set_ylabel(r"$\Delta C_\ell / C_\ell^{ref}$")
        for Cls in nn:
            l = Cls["ell"]
            ax.semilogx(l[2:], (Cls[field][2:] - ref[field][2:]) / ref[field][2:], color="r")
        ax.axhline(0, color="C0")

    return fig, ax

# def plot_Cl_rel_error(ref, Cls_nn, ax=None, field="tt", errors=None, cosmic_variance=True):
#     if not isinstance(Cls_nn, list):
#         Cls_nn = [Cls_nn]

#     if ax is None:
#         fig, ax = plt.subplots()
#     else:
#         fig = ax.figure
        

#     ax.set_xlabel(r"$\ell$")
#     ax.set_ylabel(r"$\Delta C_\ell^{{{0:}}} /  C_\ell^{{{0:}}}$".format(field.upper()))
#     ax.axhline(0, color="r", label="CLASS")

#     if errors:
#         ax.fill_between(ref["ell"][2:], ref[field][2:] - errors[2:], ref[field][2:] + errors[2:], facecolor="k", alpha=0.2)

#     # relative errors [i.e. (Cl_net - Cl_class) / Cl_class]
#     for entry in Cls_nn:
#         rel_err = (entry[field] - ref[field])[2:] / ref[field][2:]
#         ax.semilogx(ref["ell"][2:], rel_err, color="b", alpha=0.4)

#     ax.legend()

#     return fig, ax

def plot_mPk(k, mPk_ref, mPk_nn, residuals=False):
    fig, axes = plt.subplots(1 + int(residuals), sharex=True)

    if residuals:
        ax = axes[0]


    ax.set_ylabel("$P_m(k)$")
    if not residuals:
        ax.set_xlabel("$k$")

    ax.loglog(k, mPk_ref, label="CLASS full")
    ax.loglog(k, mPk_nn, label="CLASSnet", color="r", ls="--")

    ax.legend()

    if residuals:
        ax = axes[1]
        ax.semilogx(k, (mPk_nn - mPk_ref) / mPk_ref, color="r", lw=0.5)
        ax.axhline(0, color="C0")
        ax.grid()
        ax.set_xlabel("$k$")
        ax.set_ylabel(r"$\Delta P_m^(k) / P_m^{CLASS}(k)$")
        # ax.legend()
        ax.set_xlabel("$k$")

    return fig, axes

if __name__ == "__main__":
    import classy
    cosmo = classy.Class()
    cosmo.set(output="tCl,pCl,mPk")
    cosmo.compute()

    Cls = cosmo.raw_cl(2500)
    print(Cls.keys())
    import copy
    N_mock = 5
    Cls_mock = [copy.deepcopy(Cls) for _ in range(N_mock)]
    for dct in Cls_mock:
        for key in dct:
            if key == "ell":
                continue
            else:
                dct[key] *= np.random.uniform(0.9, 1.1, size=len(dct[key]))
    
    fig, axes = plt.subplots(2, sharex=True)
    plot_Cls(Cls, Cls, ax=axes[0], field="tt")
    plot_Cl_rel_error(Cls, Cls, ax=axes[1], field="tt")
    plt.show()

    import h5py as h5
    import os
    k = h5.File(
            os.path.expandvars("$CLASSNET_DATA/k_array.h5"),
            "r")["k"][()]

    mPk_ref = np.vectorize(cosmo.pk)(k, 0)

    fig, ax = plt.subplots()
    plot_mPk(k, mPk_ref, mPk_ref, ax=ax)
    plt.show()
