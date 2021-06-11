import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
import matplotlib.pyplot as plt

# style for error plots
LINESTYLE_RED = dict(ls="-", lw=0.5, alpha=0.4, color="r")
LINESTYLE_GREEN = dict(ls="-", lw=0.5, alpha=0.4, color="g")
LINESTYLE_BLUE = dict(ls="-", lw=0.5, alpha=0.4, color="b")

class SpectraPlotter:

    def __init__(self, workspace):
        self.workspace = workspace

    def plot(self, include_params=False,prefix=None,suffix=None,ylim=None):
        print("plotting spectra...")
        self._init(ylim)
        stats = self.workspace.loader().stats()
        print("building plots...")

        # this is terribly slow but whatever...
        average=[]
        for row in stats:
            av,xscale,additional = self._update_plots(row)
            average.append(av)
            last_row=row
        self._add_average(average,last_row,xscale,additional)

        if include_params:
            print("creating color plot for TT")
            self._plot_colored(stats, field="tt")
            print("creating color plot for P(k)")
            self._plot_colored(stats, field="pk")

        self._save(prefix=prefix,suffix=suffix, ylim=ylim)
        self._close_figs()

    def _init(self,ylim):
        self.figs = {
            "tt_rel_raw": plt.subplots(),
            "tt_rel_logx_raw": plt.subplots(),
            "tt_cv_raw": plt.subplots(),
            "ee_raw": plt.subplots(),
            "ee_rel_raw": plt.subplots(),
            "ee_rel_logx_raw": plt.subplots(),
            "ee_cv_raw": plt.subplots(),
            "ee_abs_raw": plt.subplots(),
            "te_raw": plt.subplots(),
            "te_lin_raw": plt.subplots(),
            "te_clean_logx_raw": plt.subplots(),
            "te_clean_lin_raw": plt.subplots(),
            "te_logx_raw": plt.subplots(),
            "te_cv_raw": plt.subplots(),
            "te_abs_raw": plt.subplots(),
            "pk": plt.subplots(),
            "pk_abs": plt.subplots(),
            "delta_m": plt.subplots(),
            "delta_m_abs": plt.subplots(),
            "pp_rel_raw": plt.subplots(),
            "pp_rel_logx_raw": plt.subplots(),
            "pp_rel_lens": plt.subplots(),
            "pp_rel_logx_lens": plt.subplots(),
            "tt_rel_lens": plt.subplots(),
            "tt_rel_logx_lens": plt.subplots(),
            "ee_lens": plt.subplots(),
            "ee_rel_lens": plt.subplots(),
            "ee_rel_logx_lens": plt.subplots(),
            "te_lens": plt.subplots(),
            "te_lin_lens": plt.subplots(),
            "te_logx_lens": plt.subplots(),
            "te_clean_lin_lens": plt.subplots(),
            "te_clean_logx_lens": plt.subplots(),
            "tt_cv_lens": plt.subplots(),
            "te_cv_lens": plt.subplots(),
            r"ee_cv_lens": plt.subplots(),
        }
        def _in_str(key,code):
            ret =False
            for i in range(len(key)-len(code)+1):
                if key[i:i+len(code)] == code:
                    ret = True
            return ret

        for _, ax in self.figs.values():
            ax.grid()

        for q in ("tt_rel_raw", "tt_rel_logx_raw", "tt_cv_raw", "te_raw","te_lin_raw" ,"te_logx_raw","te_clean_lin_raw","te_clean_logx_raw","ee_raw", "ee_rel_raw", "ee_rel_logx_raw","ee_cv_raw","te_cv_raw","pp_rel_raw","pp_rel_logx_raw",# "ee_abs_raw", "te_abs_raw",
                "tt_rel_lens", "tt_rel_logx_lens","tt_cv_lens", "te_lens","te_lin_lens","te_logx_lens","te_clean_lin_lens","te_clean_logx_lens","ee_lens", "ee_rel_lens","ee_rel_logx_lens","ee_cv_lens","te_cv_lens","pp_rel_lens","pp_rel_logx_lens"):
            fig, ax = self.figs[q]
            fig.tight_layout()
            if _in_str(q,"raw"):
                ax.set_xlim(1.5,3600)
            elif _in_str(q,"lens"):
                ax.set_xlim(1.5,3100)
            ax.set_xlabel(r"$\ell$")

        self.figs["tt_rel_raw"][1].set_ylabel(r"$\mathrm{raw} \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        self.figs["tt_rel_logx_raw"][1].set_ylabel(r"$\mathrm{raw} \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        self.figs["tt_cv_raw"][1].set_ylabel(r"$\mathrm{raw} \sqrt{\frac{2\mathrm{min}(\ell, 2000)+1}{2}} \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        ll1 = r"\ell (\ell + 1) "
        self.figs["te_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_lin_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_logx_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_clean_lin_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_clean_logx_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["ee_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"\Delta C_\ell^{EE}$")
        self.figs["ee_rel_raw"][1].set_ylabel(r"$\mathrm{raw } \Delta C_\ell^{EE}/C_\ell^{EE,true}$")
        self.figs["ee_rel_logx_raw"][1].set_ylabel(r"$\mathrm{raw } \Delta C_\ell^{EE}/C_\ell^{EE,true}$")
        self.figs["ee_raw"][1].set_yscale("symlog", linthresh=1e-16)
        self.figs["te_raw"][1].set_yscale("symlog", linthresh=1e-14)
        self.figs["ee_cv_raw"][1].set_ylabel(r"$\mathrm{raw } \left \|\Delta C_\ell^{EE} \right \|$ / (cosmic variance)")
        self.figs["te_cv_raw"][1].set_ylabel(r"$\mathrm{raw } \left \|\Delta C_\ell^{TE} \right \|$ / (cosmic variance)")
        self.figs["te_abs_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"C_\ell^{TE}$")
        self.figs["ee_abs_raw"][1].set_ylabel(r"$\mathrm{raw } "+ll1 + r"C_\ell^{EE}$")
        self.figs["ee_abs_raw"][1].set_yscale("log")

        self.figs["pk"][0].tight_layout()
        self.figs["pk"][1].set(
            xlabel=r"$k$ [Mpc${}^{-1}$]",
            ylabel=r"$(P_{NN}(k)-P_{CLASS}(k))/P_{CLASS}(k)$",
            # xlim=(1e-4, 100),
        )

        self.figs["pk_abs"][0].tight_layout()
        self.figs["pk_abs"][1].set(
            xlabel=r"$k$ [Mpc${}^{-1}$]",
            ylabel=r"$P_{NN}(k)$",
        )

        self.figs["pp_rel_raw"][1].set_ylabel(r"$\mathrm{raw} \Delta C_\ell^{PP} / C_\ell^{PP, \mathrm{true}}$")
        self.figs["pp_rel_logx_raw"][1].set_ylabel(r"$\mathrm{raw} \Delta C_\ell^{PP} / C_\ell^{PP, \mathrm{true}}$")

        self.figs["ee_lens"][1].set_yscale("symlog", linthresh=1e-16)
        self.figs["te_lens"][1].set_yscale("symlog", linthresh=1e-14)
        self.figs["te_lens"][1].set_ylabel(r"$\mathrm{lens } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_lin_lens"][1].set_ylabel(r"$\mathrm{lens } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_logx_lens"][1].set_ylabel(r"$\mathrm{lens } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_clean_lin_lens"][1].set_ylabel(r"$\mathrm{lens } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["te_clean_logx_lens"][1].set_ylabel(r"$\mathrm{lens } "+ll1 + r"\Delta C_\ell^{TE}$")
        self.figs["pp_rel_lens"][1].set_ylabel(r"$\mathrm{lensed } \Delta C_\ell^{PP} / C_\ell^{PP, \mathrm{true}}$")
        self.figs["pp_rel_logx_lens"][1].set_ylabel(r"$\mathrm{lensed } \Delta C_\ell^{PP} / C_\ell^{PP, \mathrm{true}}$")
        self.figs["tt_rel_lens"][1].set_ylabel(r"$\mathrm{lensed } \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        self.figs["tt_rel_logx_lens"][1].set_ylabel(r"$\mathrm{lensed } \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        self.figs["ee_lens"][1].set_ylabel(r"$\mathrm{lensed } "+ll1 + r"\Delta C_\ell^{EE}$")
        self.figs["ee_rel_lens"][1].set_ylabel(r"$\mathrm{lensed } \Delta C_\ell^{EE}/C_\ell^{EE,true}$")
        self.figs["ee_rel_logx_lens"][1].set_ylabel(r"$\mathrm{lensed } \Delta C_\ell^{EE}/C_\ell^{EE,true}$")
        self.figs["te_lens"][1].set_ylabel(r"$\mathrm{lensed } "+ll1 + r"\Delta C_\ell^{EE}$")
        self.figs["tt_cv_lens"][1].set_ylabel(r"$\mathrm{lensed } \sqrt{\frac{2\mathrm{min}(\ell, 2000)+1}{2}} \Delta C_\ell^{TT} / C_\ell^{TT, \mathrm{true}}$")
        self.figs["ee_cv_lens"][1].set_ylabel(r"$\mathrm{lensed } \left \|\Delta C_\ell^{EE} \right \|$ / (cosmic variance)")
        self.figs["te_cv_lens"][1].set_ylabel(r"$\mathrm{lensed } \left \|\Delta C_\ell^{TE} \right \|$ / (cosmic variance)")

        if ylim != None:
            print("\nsetting ylim\n")
            self.figs["ee_rel_raw"][1].set_ylim(-ylim["ee"],ylim["ee"])
            self.figs["ee_rel_lens"][1].set_ylim(-ylim["ee"],ylim["ee"])
            self.figs["ee_rel_logx_raw"][1].set_ylim(-ylim["ee"],ylim["ee"])
            self.figs["ee_rel_logx_lens"][1].set_ylim(-ylim["ee"],ylim["ee"])
            self.figs["tt_rel_raw"][1].set_ylim(-ylim["tt"],ylim["tt"])
            self.figs["tt_rel_lens"][1].set_ylim(-ylim["tt"],ylim["tt"])
            self.figs["tt_rel_logx_raw"][1].set_ylim(-ylim["tt"],ylim["tt"])
            self.figs["tt_rel_logx_lens"][1].set_ylim(-ylim["tt"],ylim["tt"])
            self.figs["pp_rel_raw"][1].set_ylim(-ylim["pp"],ylim["pp"])
            self.figs["pp_rel_lens"][1].set_ylim(-ylim["pp"],ylim["pp"])
            self.figs["pp_rel_logx_raw"][1].set_ylim(-ylim["pp"],ylim["pp"])
            self.figs["pp_rel_logx_lens"][1].set_ylim(-ylim["pp"],ylim["pp"])
            self.figs["te_lin_raw"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_lin_lens"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_logx_raw"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_logx_lens"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_clean_lin_raw"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_clean_lin_lens"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_clean_logx_raw"][1].set_ylim(-ylim["te"],ylim["te"])
            self.figs["te_clean_logx_lens"][1].set_ylim(-ylim["te"],ylim["te"])

    def _add_average(self, av_old,row,xscale,additional):

        cl_true = row["cl_true"]
        cl_lens_true = row["cl_lens_true"]
        ell = cl_true["ell"]
        ell_lens = cl_lens_true["ell"]
        def _in_str(key,code):
            ret =False
            for i in range(len(key)-len(code)+1):
                if key[i:i+len(code)] == code:
                    ret = True
            return ret

        l={} 
        for key in av_old[0]:
            if _in_str(key,"lens"):
                l[key]=ell_lens
            elif _in_str(key,"raw"):
                l[key]=ell
            else:
                raise ValueError("neither lensed nor raw")

        av = {}
        for key in av_old[0]:
            av[key]=np.array(l[key]-l[key],dtype="float")
        for i in range(len(av_old)):
            for key in av_old[i]:
                av[key] = av[key] + av_old[i][key]
        for key in av:
            av[key] = av[key]/len(av_old)

        for key in av:
            if key in additional:
                color = additional[key]["color"]
                lw = additional[key]["lw"]
                alpha = additional[key]["alpha"]
                key2=additional[key]["fig"]
            else:
                color="black"
                lw = 1
                alpha=1
                key2=key
            if key2 in xscale:
                if xscale[key2] == "lin":
                    self.figs[key2][1].plot(l[key2],av[key], color=color, linewidth=lw,alpha=alpha)
                else:
                    print("\nUnknown xscale key {} {}\n".format(key2, xscale[key2]))
            else:
                self.figs[key2][1].semilogx(l[key2],av[key],color=color, linewidth=lw,alpha=alpha)


    def _get_param_bounds(self, data):
        param_names = data[0]["parameters"].keys()
        params = {key: np.array([item["parameters"][key] for item in data]) for key in param_names}

        return {key: (np.min(val), np.max(val)) for key, val in params.items()}

    def _plot_colored(self, data, field="tt"):
        def plot_item_tt(ax, item, **kwargs):
            ax.set_xlabel("$\\ell$")
            ax.set_ylabel("$\\Delta C_\\ell^\mathrm{TT} / C_\\ell^\mathrm{TT}$")
            ax.semilogx(item["cl_true"]["ell"], item["cl_error_relative"]["tt"], **kwargs)

        def plot_item_pk(ax, item, **kwargs):
            ax.set_xlabel("$k$ [Mpc${}^{-1}]$")
            ax.set_ylabel("$\\Delta P(k) / P(k)$")
            ax.semilogx(item["k_pk"], item["pk_error_relative"], **kwargs)
            ax.set_yscale("symlog", linthresh=0.01)

        handlers = {
            "tt": plot_item_tt,
            "pk": plot_item_pk,
        }
        if field not in handlers:
            raise ValueError("invalid field: `{}`".format(tt))
        else:
            handler = handlers[field]

        pbounds = self._get_param_bounds(data)
        cmap = plt.cm.seismic
        for param_name, bound in pbounds.items():
            print("_plot_cl_colored({}) for {}".format(field, param_name))
            fig, ax = plt.subplots()
            norm = matplotlib.colors.Normalize(vmin=bound[0], vmax=bound[1])
            for item in data:
                param_val = item["parameters"][param_name]
                color = cmap(norm(param_val))
                handler(ax, item, color=color, lw=0.7)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            plt.colorbar(sm, label=param_name)

            ax.grid()
            dir_path = self.workspace.plots / "{}_by_param".format(field)
            dir_path.mkdir(exist_ok=True)
            path = dir_path / "{}.png".format(param_name)
            print("saving ({}, {}) to {}".format(field, param_name, path))
            fig.savefig(path, bbox_inches="tight", dpi=250)
            plt.close(fig)

    def _update_plots(self, row):
        av = {}
        xscale={}
        additional={}
        cl_true   = row["cl_true"]
        cl_nn     = row["cl_nn"]
        k_pk_true = row["k_pk"]
        pk_true   = row["pk_true"]
        k_pk_nn   = row["k_pk_nn"]
        pk_nn     = row["pk_nn"]
        ell = cl_true["ell"]
        cl_lens_true = row["cl_lens_true"]
        cl_lens_nn = row["cl_lens_nn"]
        ell_lens = cl_lens_true["ell"]

        for _, ax in self.figs.values():
            ax.axhline(0, color="k")
        
        # TT
        tt_relerr = (cl_nn["tt"] - cl_true["tt"]) / cl_true["tt"]
        self.figs["tt_rel_raw"][1].plot(ell, tt_relerr, **LINESTYLE_RED)
        av["tt_rel_raw"]=tt_relerr
        xscale["tt_rel_raw"]="lin"

        tt_lens_relerr = (cl_lens_nn["tt"] - cl_lens_true["tt"]) / cl_lens_true["tt"]
        self.figs["tt_rel_lens"][1].plot(ell_lens, tt_lens_relerr, **LINESTYLE_RED)
        av["tt_rel_lens"]=tt_lens_relerr
        xscale["tt_rel_lens"]="lin"


        # TT logx
        tt_relerr = (cl_nn["tt"] - cl_true["tt"]) / cl_true["tt"]
        self.figs["tt_rel_logx_raw"][1].semilogx(ell, tt_relerr, **LINESTYLE_RED)
        av["tt_rel_logx_raw"]=tt_relerr

        tt_lens_relerr = (cl_lens_nn["tt"] - cl_lens_true["tt"]) / cl_lens_true["tt"]
        self.figs["tt_rel_logx_lens"][1].semilogx(ell_lens, tt_lens_relerr, **LINESTYLE_RED)
        av["tt_rel_logx_lens"]=tt_lens_relerr


        cosmic_variance = np.sqrt(2 / (2 * np.minimum(ell, 2000) + 1))
        tt_relerr_cv = tt_relerr / cosmic_variance
        self.figs["tt_cv_raw"][1].semilogx(ell, tt_relerr_cv, **LINESTYLE_RED)
        av["tt_cv_raw"]=tt_relerr_cv

        cosmic_variance_lens = np.sqrt(2 / (2 * np.minimum(ell_lens, 2000) + 1))
        tt_lens_relerr_cv = tt_lens_relerr / cosmic_variance_lens
        self.figs["tt_cv_lens"][1].semilogx(ell_lens, tt_lens_relerr_cv, **LINESTYLE_RED)
        av["tt_cv_lens"]=tt_lens_relerr_cv

        def plot_err(ax, qty, style):
            ax.semilogx(ell, ell * (ell + 1) * qty, **style)

        def plot_err_pm(ax, qty, style):
            plot_err(ax, qty, style)
            plot_err(ax, -qty, style)
            
        def plot_err_lens(ax, qty, style):
            ax.semilogx(ell_lens, ell_lens * (ell_lens + 1) * qty, **style)

        def plot_err_pm_lens(ax, qty, style):
            plot_err_lens(ax, qty, style)
            plot_err_lens(ax, -qty, style)
        
        pp_relerr = (cl_nn["pp"] - cl_true["pp"]) / cl_true["pp"]
        self.figs["pp_rel_raw"][1].plot(ell, pp_relerr, **LINESTYLE_RED)
        av["pp_rel_raw"]=pp_relerr
        xscale["pp_rel_raw"]="lin"
        pp_lens_relerr = (cl_lens_nn["pp"] - cl_lens_true["pp"]) / cl_lens_true["pp"]
        self.figs["pp_rel_lens"][1].plot(ell_lens, pp_lens_relerr, **LINESTYLE_RED)
        av["pp_rel_lens"]=pp_lens_relerr
        xscale["pp_rel_lens"]="lin"

        pp_relerr = (cl_nn["pp"] - cl_true["pp"]) / cl_true["pp"]
        self.figs["pp_rel_logx_raw"][1].semilogx(ell, pp_relerr, **LINESTYLE_RED)
        av["pp_rel_logx_raw"]=pp_relerr
        pp_lens_relerr = (cl_lens_nn["pp"] - cl_lens_true["pp"]) / cl_lens_true["pp"]
        self.figs["pp_rel_logx_lens"][1].semilogx(ell_lens, pp_lens_relerr, **LINESTYLE_RED)
        av["pp_rel_logx_lens"]=pp_lens_relerr

        # TE + EE
        for q in ("ee", "te"):
            err = (cl_nn[q] - cl_true[q])
            # err_relmax = err / cl_true[q].max()
            ax = self.figs[q+"_raw"][1]

            # cosmic variance
            cv = np.sqrt(2 / (2 * ell + 1)) * cl_true[q]
            ll1_cv = ell * (ell + 1) * cv
            ax.semilogx(ell, ll1_cv, lw=0.4, color="purple")
            ax.semilogx(ell, -ll1_cv, lw=0.4, color="purple")

            # 1% and 0.1%
            ls_blue = LINESTYLE_BLUE.copy()
            ls_blue["alpha"] = 0.4
            ls_green = LINESTYLE_GREEN.copy()
            ls_green["alpha"] = 0.4
            plot_err_pm(ax, cl_true[q] /  100.0, ls_blue)
            plot_err_pm(ax, cl_true[q] / 1000.0, ls_green)

            # actual error
            plot_err(ax, err, LINESTYLE_RED)
            av[q+"_raw"]=ell*(ell+1)*err

            ###### relative to cosmic variance ####
            ax = self.figs[q + "_cv_raw"][1]
            ax.loglog(ell, np.abs(err / cv), **LINESTYLE_RED)
            av[q+"_cv_raw"]=np.abs(err/cv)

        # ee_rel_raw
        ax = self.figs["ee_rel_raw"][1]
        err = cl_nn["ee"] - cl_true["ee"]
        ax.plot(ell,err/cl_true["ee"],**LINESTYLE_RED)
        av["ee_rel_raw"]=err/cl_true["ee"]
        xscale["ee_rel_raw"]="lin"
        # ee_rel_lens
        ax = self.figs["ee_rel_lens"][1]
        err = cl_lens_nn["ee"] - cl_lens_true["ee"]
        ax.plot(ell_lens,err/cl_lens_true["ee"],**LINESTYLE_RED)
        av["ee_rel_lens"]=err/cl_lens_true["ee"]
        xscale["ee_rel_lens"]="lin"
        # ee_rel_logx_raw
        ax = self.figs["ee_rel_logx_raw"][1]
        err = cl_nn["ee"] - cl_true["ee"]
        ax.semilogx(ell,err/cl_true["ee"],**LINESTYLE_RED)
        av["ee_rel_logx_raw"]=err/cl_true["ee"]
        # ee_rel_logx_lens
        ax = self.figs["ee_rel_logx_lens"][1]
        err = cl_lens_nn["ee"] - cl_lens_true["ee"]
        ax.semilogx(ell_lens,err/cl_lens_true["ee"],**LINESTYLE_RED)
        av["ee_rel_logx_lens"]=err/cl_lens_true["ee"]
        
        
        # TE rel_raw
        err = (cl_nn["te"] - cl_true["te"])
        ax = self.figs["te_lin_raw"][1]
        # 1% and 0.1%
        av["tem_10_lin_raw"]=-np.abs(ell * (ell+ 1) * cl_true["te"]/10.0)
        additional["tem_10_lin_raw"]={"fig":"te_lin_raw","lw":0.8, "alpha":0.2, "color":"purple"}
        av["tem_100_lin_raw"]=-np.abs(ell * (ell+ 1) * cl_true["te"]/100.0)
        additional["tem_100_lin_raw"]={"fig":"te_lin_raw","lw":0.8, "alpha":0.2, "color":"blue"}
        av["tem_1000_lin_raw"]=-np.abs(ell * (ell+ 1) * cl_true["te"]/1000.0)
        additional["tem_1000_lin_raw"]={"fig":"te_lin_raw","lw":0.8, "alpha":0.2, "color":"green"}
        av["te_10_lin_raw"]=np.abs(ell * (ell+ 1) * cl_true["te"]/10.0)
        additional["te_10_lin_raw"]={"fig":"te_lin_raw","lw":0.8, "alpha":0.2, "color":"purple"}
        av["te_100_lin_raw"]=np.abs(ell * (ell+ 1) * cl_true["te"]/100.0)
        additional["te_100_lin_raw"]={"fig":"te_lin_raw","lw":0.8, "alpha":0.2, "color":"blue"}
        av["te_1000_lin_raw"]=np.abs(ell * (ell+ 1) * cl_true["te"]/1000.0)
        additional["te_1000_lin_raw"]={"fig":"te_lin_raw","lw":0.8, "alpha":0.2, "color":"green"}
        # actual error
        ax.plot(ell,ell*(ell+1)*err, **LINESTYLE_RED)
        av["te_lin_raw"]=ell*(ell+1)*err
        xscale["te_lin_raw"]="lin"
   
        # TE rel_lens
        err_lens = (cl_lens_nn["te"] - cl_lens_true["te"])
        ax = self.figs["te_lin_lens"][1]
        # 1% and 0.1%
        av["tem_10_lin_lens"]=-np.abs(ell_lens * (ell_lens+ 1) * cl_lens_true["te"]/10.0)
        additional["tem_10_lin_lens"]={"fig":"te_lin_lens","lw":0.8, "alpha":0.2, "color":"purple"}
        av["tem_100_lin_lens"]=-np.abs(ell_lens * (ell_lens+ 1) * cl_lens_true["te"]/100.0)
        additional["tem_100_lin_lens"]={"fig":"te_lin_lens","lw":0.8, "alpha":0.2, "color":"blue"}
        av["tem_1000_lin_lens"]=-np.abs(ell_lens * (ell_lens+ 1) * cl_lens_true["te"]/1000.0)
        additional["tem_1000_lin_lens"]={"fig":"te_lin_lens","lw":0.8, "alpha":0.2, "color":"green"}
        av["te_10_lin_lens"]=np.abs(ell_lens * (ell_lens+ 1) * cl_lens_true["te"]/10.0)
        additional["te_10_lin_lens"]={"fig":"te_lin_lens","lw":0.8, "alpha":0.2, "color":"purple"}
        av["te_100_lin_lens"]=np.abs(ell_lens * (ell_lens+ 1) * cl_lens_true["te"]/100.0)
        additional["te_100_lin_lens"]={"fig":"te_lin_lens","lw":0.8, "alpha":0.2, "color":"blue"}
        av["te_1000_lin_lens"]=np.abs(ell_lens * (ell_lens+ 1) * cl_lens_true["te"]/1000.0)
        additional["te_1000_lin_lens"]={"fig":"te_lin_lens","lw":0.8, "alpha":0.2, "color":"green"}
        # actual error
        ax.plot(ell_lens,ell_lens*(ell_lens+1)*err_lens, **LINESTYLE_RED)
        av["te_lin_lens"]=ell_lens*(ell_lens+1)*err_lens
        xscale["te_lin_lens"]="lin"
        
        # TE rel_raw
        err = (cl_nn["te"] - cl_true["te"])
        ax = self.figs["te_logx_raw"][1]
        # 1% and
        av["tem_10_logx_raw"]=-np.abs(ell * (ell+ 1) * cl_true["te"]/10.0)
        additional["tem_10_logx_raw"]={"fig":"te_logx_raw","lw":0.8, "alpha":0.2, "color":"purple"}
        av["tem_100_logx_raw"]=-np.abs(ell * (ell+ 1) * cl_true["te"]/100.0)
        additional["tem_100_logx_raw"]={"fig":"te_logx_raw","lw":0.8, "alpha":0.2, "color":"blue"}
        av["tem_1000_logx_raw"]=-np.abs(ell * (ell+ 1) * cl_true["te"]/1000.0)
        additional["tem_1000_logx_raw"]={"fig":"te_logx_raw","lw":0.8, "alpha":0.2, "color":"green"}
        av["te_10_logx_raw"]=np.abs(ell * (ell+ 1) * cl_true["te"]/10.0)
        additional["te_10_logx_raw"]={"fig":"te_logx_raw","lw":0.8, "alpha":0.2, "color":"purple"}
        av["te_100_logx_raw"]=np.abs(ell * (ell+ 1) * cl_true["te"]/100.0)
        additional["te_100_logx_raw"]={"fig":"te_logx_raw","lw":0.8, "alpha":0.2, "color":"blue"}
        av["te_1000_logx_raw"]=np.abs(ell * (ell+ 1) * cl_true["te"]/1000.0)
        additional["te_1000_logx_raw"]={"fig":"te_logx_raw","lw":0.8, "alpha":0.2, "color":"green"}
        # actual error
        ax.semilogx(ell,ell*(ell+1)*err, **LINESTYLE_RED)
        av["te_logx_raw"]=ell*(ell+1)*err
        
        # TE rel_lens
        err_lens = (cl_lens_nn["te"] - cl_lens_true["te"])
        ax = self.figs["te_logx_lens"][1]
        # 1% and 0.1%
        av["tem_10_logx_lens"]=-np.max(ell_lens * (ell_lens + 1) * cl_lens_true["te"]/10.0)
        additional["tem_10_logx_lens"]={"fig":"te_logx_lens","lw":0.8, "alpha":0.2, "color":"purple"}
        av["tem_100_logx_lens"]=-np.max(ell_lens * (ell_lens + 1) * cl_lens_true["te"]/100.0)
        additional["tem_100_logx_lens"]={"fig":"te_logx_lens","lw":0.8, "alpha":0.2, "color":"blue"}
        av["tem_1000_logx_lens"]=-np.max(ell_lens * (ell_lens + 1) * cl_lens_true["te"]/1000.0)
        additional["tem_1000_logx_lens"]={"fig":"te_logx_lens","lw":0.8, "alpha":0.2, "color":"green"}
        av["te_10_logx_lens"]=np.max(ell_lens * (ell_lens + 1) * cl_lens_true["te"]/10.0)
        additional["te_10_logx_lens"]={"fig":"te_logx_lens","lw":0.8, "alpha":0.2, "color":"purple"}
        av["te_100_logx_lens"]=np.max(ell_lens * (ell_lens + 1) * cl_lens_true["te"]/100.0)
        additional["te_100_logx_lens"]={"fig":"te_logx_lens","lw":0.8, "alpha":0.2, "color":"blue"}
        av["te_1000_logx_lens"]=np.max(ell_lens * (ell_lens + 1) * cl_lens_true["te"]/1000.0)
        additional["te_1000_logx_lens"]={"fig":"te_logx_lens","lw":0.8, "alpha":0.2, "color":"green"}
        # actual error
        ax.semilogx(ell_lens,ell_lens*(ell_lens+1)*err_lens, **LINESTYLE_RED)
        av["te_logx_lens"]=ell_lens*(ell_lens+1)*err_lens
        
        
        err = (cl_nn["te"] - cl_true["te"])
        ax = self.figs["te_clean_lin_raw"][1]
        ax.plot(ell,ell*(ell+1)*err, **LINESTYLE_RED)
        av["te_clean_lin_raw"]=ell*(ell+1)*err
        xscale["te_clean_lin_raw"]="lin"

        err_lens = (cl_lens_nn["te"] - cl_lens_true["te"])
        ax = self.figs["te_clean_lin_lens"][1]
        ax.plot(ell_lens,ell_lens*(ell_lens+1)*err_lens, **LINESTYLE_RED)
        av["te_clean_lin_lens"]=ell_lens*(ell_lens+1)*err_lens
        xscale["te_clean_lin_lens"]="lin"

        err = (cl_nn["te"] - cl_true["te"])
        ax = self.figs["te_clean_logx_raw"][1]
        ax.semilogx(ell,ell*(ell+1)*err, **LINESTYLE_RED)
        av["te_clean_logx_raw"]=ell*(ell+1)*err
        
        err_lens = (cl_lens_nn["te"] - cl_lens_true["te"])
        ax = self.figs["te_clean_logx_lens"][1] 
        ax.semilogx(ell_lens,ell_lens*(ell_lens+1)*err_lens, **LINESTYLE_RED)
        av["te_clean_logx_lens"]=ell_lens*(ell_lens+1)*err_lens


        for q in ("ee", "te"):
            err_lens = (cl_lens_nn[q] - cl_lens_true[q])
            # err_relmax = err / cl_true[q].max()
            ax = self.figs[q+"_lens"][1]

            # cosmic variance
            cv_lens = np.sqrt(2 / (2 * ell_lens + 1)) * cl_lens_true[q]
            ll1_cv_lens = ell_lens * (ell_lens + 1) * cv_lens
            ax.semilogx(ell_lens, ll1_cv_lens, lw=0.4, color="purple")
            ax.semilogx(ell_lens, -ll1_cv_lens, lw=0.4, color="purple")

            # 1% and 0.1%
            ls_blue = LINESTYLE_BLUE.copy()
            ls_blue["alpha"] = 0.4
            ls_green = LINESTYLE_GREEN.copy()
            ls_green["alpha"] = 0.4
            plot_err_pm_lens(ax, cl_lens_true[q] /  100.0, ls_blue)
            plot_err_pm_lens(ax, cl_lens_true[q] / 1000.0, ls_green)

            # actual error
            plot_err_lens(ax, err_lens, LINESTYLE_RED)
            av[q+"_lens"]=ell_lens*(ell_lens+1)*err_lens


            ###### relative to cosmic variance ####
            ax = self.figs[q + "_cv_lens"][1]
            ax.loglog(ell_lens, np.abs(err_lens / cv_lens), **LINESTYLE_RED)
            av[q+"_cv_lens"]=np.abs(err_lens/cv_lens)

        # TE + EE absolute
        for q in ("ee", "te"):
            ax = self.figs[q + "_abs_raw"][1]
            ax.semilogx(ell, ell * (ell + 1) * cl_true[q], **LINESTYLE_GREEN)
            ax.semilogx(ell, ell * (ell + 1) * cl_nn[q], **LINESTYLE_RED)

        # P(k)
        # since P_NN(k) and P_true(k) may be sampled on different k grids, we
        # need to interpolate (in this case, onto the k_pk_true)
        REINTERP_PK = True
        if REINTERP_PK:
            pk_spline = CubicSpline(k_pk_nn, pk_nn)
            pk_nn_resampled = pk_spline(k_pk_true)
            pk_relerr = (pk_nn_resampled - pk_true) / pk_true
        else:
            assert np.allclose(k_pk_nn, k_pk_true)
            pk_relerr = (pk_nn - pk_true) / pk_true
        self.figs["pk"][1].semilogx(k_pk_true, pk_relerr, **LINESTYLE_RED)
        self.figs["pk"][1].set_yscale("symlog", linthresh=0.001)

        # this will raise a warning that there are no positive values and
        # hence loglog is not possible in matplotlib version 3.3.1 but this
        # appears to be a bug in matplotlib and therefore can be ignored.
        self.figs["pk_abs"][1].loglog(k_pk_nn, pk_nn, **LINESTYLE_RED)
        self.figs["pk_abs"][1].loglog(k_pk_true, pk_true, **LINESTYLE_GREEN)

        # delta_m
        self.figs["delta_m_abs"][1].loglog(row["k"], -row["delta_m"], color="g")
        self.figs["delta_m_abs"][1].loglog(row["k_nn"], -row["delta_m_nn"], color="r")

        dm_nn_interp = CubicSpline(row["k_nn"], row["delta_m_nn"])(row["k"])
        dm_rel_err = (dm_nn_interp - row["delta_m"]) / row["delta_m"]
        self.figs["delta_m"][1].semilogx(row["k"], dm_rel_err, **LINESTYLE_RED)
        self.figs["delta_m"][1].set_yscale("symlog", linthresh=0.001)
        
        return av,xscale, additional

    def _save(self, prefix=None,suffix=None,ylim=None):
        if suffix!=None and suffix!="":
            suffix="_"+suffix
        for name, (fig, _) in self.figs.items():
            if ylim != None and name[:2] in ["ee","tt","te","pp"]:
                ylimits = "_ylim_"+str(ylim[name[:2]])
                print("save_ylim")
            else:
                ylimits = ""
            if not prefix:
                path = self.workspace.plots / "{}{}{}.pdf".format(name,suffix,ylimits)
            else:
                dir_path = self.workspace.plots / prefix
                dir_path.mkdir(parents=True, exist_ok=True)
                path = dir_path / "{}{}{}.pdf".format(name,suffix,ylimits)

            print("saving plot to", path)
            fig.savefig(path, bbox_inches="tight")

    def _close_figs(self):
        for (fig, _) in self.figs.values():
            plt.close(fig)
        self.figs = dict()
