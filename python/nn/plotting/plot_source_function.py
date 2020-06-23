import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.ticker as ticker
from scipy.interpolate import interp2d

def fmt_colorbar(x, pos):
    a, b = '{: .2e}'.format(x).split('e')
    a = a.replace(" ", "\\ ")
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def plot_source_function(ax, k, tau, S, tau_rec=None, levels=256, cmap="RdBu_r", 
        title=None, title_in_plot=True, 
        use_contourf=True,
        show_xlabel=True, show_ylabel=True,
        show_colorbar=True,
        return_cont=False,
        vmax=None):
    # Find maximum absolute value of S in order to center the center of the 
    # color scale (e.g. white) to a value of S = 0
    vmax = vmax if vmax is not None else np.max(np.abs(S))
    # Prepare tau for plotting, i.e. renormalize it w.r.t. to tau_rec, if given
    tau_plot = tau if tau_rec is None else tau / tau_rec

    # Set up axis scaling
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Perform actual plotting; either using contourf or imshow, depending on 
    # given argument `use_contourf`
    if use_contourf:
        contours = ax.contourf(k, tau_plot, S.transpose(), levels=levels, cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.invert_yaxis()
    else:
        contours = ax.imshow(S.transpose(), 
                vmin=-vmax, vmax=vmax, 
                cmap=cmap, 
                origin="upper", 
                aspect="auto",
                extent=[k.min(), k.max(), tau_plot.min(), tau_plot.max()])

    # Axis labels
    if show_xlabel:
        ax.set_xlabel("$k [{Mpc}^{-1}$]")
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()

    if show_ylabel:
        if tau_rec is None:
            ax.set_ylabel("$\\tau$ [Mpc${}^{-1}$]")
        else:
            ax.set_ylabel("$\\tau / \\tau_{rec}$")
    if title:
        if not title_in_plot:
            ax.set_title(title)
        else:
            ax.text(0.05, 0.05, title, transform=ax.transAxes)

    if show_colorbar:
        # Adjust colorbar to show exponential format
        # colorbar_formatter = ticker.FuncFormatter(fmt_colorbar)
        # cb = plt.colorbar(contours, ax=ax, format=colorbar_formatter)
        cb = plt.colorbar(contours, ax=ax)
        cb.formatter.set_powerlimits((0, 0))
        cb.ax.yaxis.set_offset_position('left')                         
        cb.update_ticks()

    if return_cont:
        return ax, contours
    else:
        return ax

def plot_source_functions(
        k, tau, sources, 
        k_nn, tau_nn, sources_nn, 
        tau_rec=None,
        size=3, 
        levels=50,
        selection=None, labels=None, 
        transpose_plot=False,
        plot_losses=False,
        losses=None,
        val_losses=None,
        epochs=None,
        loss_fields=None,
        plot_Cls=False,
        Cls=None, Cls_nn=None,
        ):
    """
    Plot two sets of source functions and compare.

    :param k: k sampling for `sources`
    :param tau: tau sampling for `sources`
    :param sources: dictionary containing source functions obtained without NNs

    :param k_nn: k sampling for `sources_nn`
    :param tau_nn: tau sampling for `sources_nn`
    :param sources_nn: dictionary containing source functions obtained with NNs

    :param selection: List of names of source functions (i.e. keys into `sources` / `sources_nn`) to plot
    :param labels: List of labels; Must contain a %s placeholder which will accept either "CLASS" or "net"
    :param size: size (in inches) of each source function plot
    :param levels: number of contour levels
    :param transpose_plot: boolean; if `False`, source functions will be plotted in columns; otherwise, in rows.

    :param plot_losses: If True, add a plot for each source function showing the losses of the corresponding model
    :param losses: dictionary mapping source function name -> array of losses
    :param val_losses: dictionary mapping source function name -> array of val_losses
    :param epochs: dictionary mapping source function name -> current epoch (used to highlight current point)

    :returns: created figure
    """

    if selection is None:
        selection = sorted(set(sources) & set(sources_nn))
    if labels is None:
        labels = selection

    if plot_Cls:
        if not Cls or not Cls_nn:
            raise ValueError("Need to provide Cls!")

    if plot_losses:
        assert losses is not None
        if loss_fields is None:
            loss_fields = list(losses)

    if transpose_plot:
        nrows = len(selection) + int(plot_Cls)
        ncols = 3 + int(plot_losses)
    else:
        nrows = 3 + int(plot_losses)
        ncols = len(selection) + int(plot_Cls)

    figsize = (ncols * (1.8 * size), nrows * size)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    
    proxy = axes.T if transpose_plot else axes


    for i, (fun, label) in enumerate(zip(selection, labels)):
        if transpose_plot:
            show_xlabel = i == 0
            show_ylabel = True
        else:
            show_xlabel = True
            show_ylabel = i == 0

        plot_source_function(
                proxy[0, i], 
                k, tau, sources[fun], 
                tau_rec=tau_rec, 
                levels=levels, 
                show_xlabel=show_xlabel, show_ylabel=show_ylabel,
                title=label % "CLASS")

        if transpose_plot:
            show_ylabel = False
        else:
            show_xlabel = False

        plot_source_function(
                proxy[1, i], 
                k_nn, tau_nn, sources_nn[fun], 
                tau_rec=tau_rec, 
                levels=levels, 
                show_xlabel=show_xlabel, show_ylabel=show_ylabel,
                title=label % "net")
        
        spline = interp2d(tau, k, sources[fun])
        source_on_nn_grid = spline(tau_nn, k_nn)
        diff = sources_nn[fun] - source_on_nn_grid

        plot_source_function(
                proxy[2, i], 
                k_nn, tau_nn, diff, 
                tau_rec=tau_rec, 
                levels=levels, 
                show_xlabel=show_xlabel, show_ylabel=show_ylabel,
                title="$\\Delta$")

        if plot_losses:
            ax = proxy[3, i]
            if fun in loss_fields:
                ax.semilogy(losses[fun], label="Loss", marker="o", color="C3")
                xtick_skip = max(int(len(losses[fun]) / 10), 1)
                xticks = np.arange(0, len(losses[fun]), xtick_skip)
                ax.set_xticks(xticks)
                if val_losses is not None:
                    ax.semilogy(val_losses[fun], label="Val. Loss", marker="o", color="C2")
                if epochs:
                    ax.axvline(epochs[fun], color="k", ls=":", lw=2, label="current epoch")
                ax.set_xlabel("Epoch")
                ax.legend()
                ax.grid()
            else:
                ax.axis("off")

    # if plot_Cls:
    #     ax = proxy[0, i]
    #     ax.set_xlabel(r"$\ell$")
    #     ax.set_ylabel(r"$\ell (\ell + 1) C_\ell^{{{}}} / (2 \pi)$".format(field.upper()))
    #     if errors:
    #         ax.fill_between(ref["ell"][2:], ref[field][2:] - errors[2:], ref[field][2:] + errors[2:], facecolor="k", alpha=0.2)

    # for ax in axes.flat:
    #     ax.label_outer()
    # fig.subplots_adjust(hspace=0, wspace=0.4)

    return fig

def animate_source_function_training(k, tau, Ss, **kwargs):
    fig, ax = plt.subplots()

    vmax = np.max(np.abs(Ss[-1]))

    _, cont = plot_source_function(ax, k, tau, Ss[0], vmax=vmax, return_cont=True, show_colorbar=False, **kwargs)
    def plot_frame(frame_idx):
        print("anim: on frame {} / {}".format(frame_idx + 1, len(Ss)))
        _, cont = plot_source_function(ax, k, tau, Ss[frame_idx], vmax=vmax, return_cont=True, show_colorbar=False, **kwargs)
        return cont,

    # anim = animation.FuncAnimation(fig, plot_frame, len(Ss))
    anim = animation.FuncAnimation(fig, plot_frame, 1)
    return anim

def plot_source_function_training_history(k, tau, Ss, size=4, **kwargs):
    fig, axes = plt.subplots(nrows=len(Ss), figsize=(size, len(Ss) * size))
    vmax = np.max(np.abs(Ss[-1]))
    for S, ax in zip(Ss, axes):
        plot_source_function(ax, k, tau, S, vmax=vmax)
    return fig

