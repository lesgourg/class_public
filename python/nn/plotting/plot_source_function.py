import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.ticker as ticker
import matplotlib
from scipy.interpolate import interp2d, CubicSpline
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('font',**{'family':'serif','serif':['cm'],'size':15})
rc('text', usetex=True)
plt.rcParams['axes.titley']=1.03
plt.rcParams["figure.autolayout"]=True
from scipy.interpolate import interp2d

SOURCE_LATEX = {
    "t0_reco_no_isw": r"T_0,\mathrm{reco}",
    "t0_reio_no_isw": r"T_0,\mathrm{reio}",
    "t0_isw": r"T_0,\mathrm{ISW}",
    "t1": r"T_1",
    "t2_reco": r"T_2,\mathrm{reco}",
    "t2_reio": r"T_2,\mathrm{reio}",
    "phi_plus_psi": r"\phi+\psi",
    "delta_m": r"\delta_m",
    "delta_cb": r"\delta_{cb}",
}

def fmt_colorbar(x, pos):
    a, b = '{: .2e}'.format(x).split('e')
    a = a.replace(" ", "\\ ")
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def plot_source_function(ax, S, k, tau, S_name=None, 
        use_nn=True,
        levels=256, 
        cmap="PuOr",
        title="auto",
        show_xlabel=True, show_ylabel=True,
        show_colorbar=True,
        return_cont=False,
        vmax=None):

    ax.loglog()

    xx, yy = np.meshgrid(tau,k)
    shading="gouraud"
    if S_name in ["phi_plus_psi","delta_cb","delta_m"]:
        if S[0][0]>0:
            minS = min([abs(np.min(S)),abs(np.min(S))])
        else:
            minS = -max([abs(np.max(S)),abs(np.max(S))])
        im = ax.pcolormesh(xx,yy,S, shading=shading, rasterized=True, 
                norm=matplotlib.colors.SymLogNorm(linthresh=abs(minS)),
                cmap=matplotlib.colors.LinearSegmentedColormap.from_list("Custom",plt.get_cmap(cmap)(np.linspace(0.25+0.25*minS/abs(minS),0.75+0.25*minS/abs(minS),50)))
                )
    else:
        max_absS = np.max([abs(np.max(S)),abs(np.min(S))])
        im = ax.pcolormesh(xx,yy,S, cmap=cmap, shading=shading, rasterized=True, vmin=-max_absS,vmax=max_absS)

    # Axis labels
    if show_xlabel:
        ax.set_xlabel("$\\tau$ [Mpc${}^{-1}$]")
    if show_ylabel:
        ax.set_ylabel("$k$ [Mpc${}^{-1}$]")
    if title=="auto":
        if use_nn:
            my_title = r"$S^\mathrm{Net}_{%(x)s}$"%{"x":SOURCE_LATEX[S_name]}
        else:
            my_title = r"$S^\mathrm{Full}_{%(x)s}$"%{"x":SOURCE_LATEX[S_name]}
        ax.set_title(my_title)
    else:
        ax.set_title(title)


    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    if show_colorbar:
        cb = plt.colorbar(im, ax=ax)



def plot_source_function_difference(axs, S_name, 
        S_NN, k_NN, tau_NN, 
        S_FULL, k_FULL, tau_FULL, 
        levels=256, 
        cmap="PuOr",
        title=True,
        show_xlabel=True, show_ylabel=True,
        show_colorbar=True,
        return_cont=False,
        vmax=None):

    # The source functions obtained by CLASSFULL have to be interpolated.
    # Due to the logarithimc character of ['phi_plus_psi','delta_cb','delta_m'],
    # we are forced to use two different interpolations
    if S_name in ['phi_plus_psi','delta_cb','delta_m']:
        S_FULL_IP = interp2d(tau_FULL, k_FULL, S_FULL)(tau_NN, k_NN)
        S_DIFF_REL = (S_NN-S_FULL_IP)/S_FULL_IP
        my_title = r"$(S^\mathrm{Net}_{%(x)s}-S^\mathrm{Full}_{%(x)s})/S^\mathrm{Full}_{%(x)s}$"%{"x":SOURCE_LATEX[S_name]}
    else:
        S_FULL_IP = np.array([CubicSpline(k_FULL,S_FULL[:,itau])(k_NN) for itau in range(len(tau_FULL))]).T
        S_DIFF_REL = S_NN-S_FULL_IP
        my_title = r"$S^\mathrm{Net}_{%(x)s}-S^\mathrm{Full}_{%(x)s}$"%{"x":SOURCE_LATEX[S_name]}

    # First plot individual contributions
    plot_source_function(axs[0],S_FULL_IP,k_NN, tau_NN, S_name=S_name, use_nn=False)
    plot_source_function(axs[1],S_NN,k_NN, tau_NN, S_name=S_name, use_nn=True)

    # Generate the title for the difference plot

    # Then plot the relative difference
    plot_source_function(axs[2],S_DIFF_REL,k_NN, tau_NN, title=my_title)


