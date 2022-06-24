import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import pkg_resources

import scipy as sc

from classynet.generate.generate_spectra import PK_M_AT_Z,PK_CB_AT_Z

# Settings for different spectra
REQUIERES_ELL_FACTOR = ['tt','bb','ee','te','tp']
REQUIERES_ELL_SQUARE_FACTOR = ['pp']
SKIP_REL_DIF = ['te','tp']
LOG_AXIS = ['pk','pk_cb','tp']
SYMLOG_AXIS = ['te']
SYMLOG_LIN_THRESHOLD = {'te':1e-12,'pp':1e-11}
LOG_XAXIS = ['ee','te','pk','pk_cb','tt','pp', 'bb']
HAS_COSMIC_VARIANCE = ['tt','pp','ee','te','bb']
HAS_NOISE = ['pp','tt','ee','te']
LOG_YLIMITS= {'tt':1e-13,'pp':1e-12,'ee':1e-16,'bb':1e-22}
HAS_ELL = ['ee','te','tt','bb']
HAS_ELL2 = ['pp']

# latex dict
LATEX_DICT = {
    'tt':r'\mathrm{TT}',
    'bb':r'\mathrm{BB}',
    'ee':r'\mathrm{EE}',
    'te':r'\mathrm{TE}',
    'pp':r'\phi \phi',
    'pk':r'P_{\mathrm{m}}',
    'pk_cb':r'P_{\mathrm{cb}}',
    'tp':r'\mathrm{T} \phi'
}

# latex dict parameters
LATEX_DICT_PARAMETER = {
    'omega_b':r'\Omega_\mathrm{b} h^2',
    'omega_cdm':r'\Omega_\mathrm{c} h^2',
    'omega_m':r'\Omega_\mathrm{m} h^2',
    'h':r'h',
    'tau_reio':r'\kappa_\mathrm{reio}',
    'w0_fld':r'w_0',
    'wa_fld':r'w_a',
    'N_ur':r'\triangle N_\mathrm{eff}',
    'omega_ncdm':r'\Omega_\mathrm{ncdm} h^2',
    'Omega_k':r'\Omega_\mathrm{k}',
}



#
# The following class is a collection of different routines.
# It will be called numerous times in the tester class.
#
class CL_plotter:

    def __init__(self,
        figsize=(7,5),
        linewidth = 0.5):
        self.figsize = figsize
        self.linewidth = linewidth

        rc('font',**{'family':'serif','serif':['cm'],'size':15})
        rc('text', usetex=True)


    def cut_datasets(self,FULL_cls,NN_cls,spectrum,N):
        """
        This function takes the input dicts and cuts 
        """
        if spectrum in ['pk','pk_cb']:
            # Load data
            x = np.array(FULL_cls['kk'])
            y_FULL = np.array(FULL_cls[spectrum])
            y_NN = np.array(NN_cls[spectrum])

            # cut where pk=0
            y_FULL=y_FULL[:N]
            y_NN=y_NN[:N]

            #trivial factor
            y_factor = np.ones(len(x[0]))
        else:
            # load data
            x = np.repeat([np.arange(len(FULL_cls[spectrum][0])-2) + 2],N,axis=0)
            y_FULL = np.array(FULL_cls[spectrum])
            y_NN = np.array(NN_cls[spectrum])

            # cut first indices
            y_FULL = y_FULL[:N,2:]
            y_NN = y_NN[:N,2:]

            if spectrum in HAS_ELL2:
                y_factor = (x[0] * (x[0] + 1))**2
            elif spectrum in HAS_ELL:
                y_factor = x[0] * (x[0] + 1)
            else:
                y_factor = np.ones(len(x[0]))

        return(x,y_FULL,y_NN,y_factor)

    def calculate_cosmic_variance(self, ell, spectrum):

        if spectrum in ["ee","tt","pp","bb"]:
            cv = np.sqrt(2/(2*ell+1))
            return cv

    def calculate_cosmic_variance_te(self, ell, te, tt, ee):
        return()

    def load_planck_noise(self,FULL_cls):
        resource_path = '/'.join(('plotting', 'approx_planck_noise.dat'))
        noise_path = pkg_resources.resource_stream('classynet', resource_path)
        noise = np.genfromtxt(noise_path)
        l_pp = noise[:1998,0]

        # calculate relevant cosmic variances
        tt_cv = self.calculate_cosmic_variance(noise[:,0],'tt')
        ee_cv = self.calculate_cosmic_variance(noise[:,0],'ee')
        pp_cv = self.calculate_cosmic_variance(l_pp,'pp')

        # calculate different noise contributions
        tt_noise = np.sqrt(2/(2*FULL_cls["ell"][0][2:]+1))*noise[:,1]/FULL_cls['tt'][0][2:]/(2.726*10**6)**2
        ee_noise = np.sqrt(2/(2*FULL_cls["ell"][0][2:]+1))*noise[:,2]/FULL_cls['ee'][0][2:]/(2.726*10**6)**2

        # transform deflection spectrum into lensing potential        
        dd_noise = noise[:1998,3] / ((l_pp+1)*l_pp)**2 * 2 * np.math.pi #/(2.726*10**6)**2

        # calculate the relative lensing noise
        pp_noise = np.sqrt(2/(2*l_pp+1))*dd_noise/FULL_cls['pp'][0][2:2000]

        # caclualte TE component
        cv_term = FULL_cls['te'][0][2:] * FULL_cls['te'][0][2:]
        noise_cv_term = (FULL_cls["tt"][0][2:] + np.array(noise[:,1])/(2.726*10**6)**2)  *  (FULL_cls["ee"][0][2:] + np.array(noise[:,2])/(2.726*10**6)**2)
        te_noise = np.sqrt(1/(2*FULL_cls["ell"][0][2:]+1)) * np.sqrt(cv_term + noise_cv_term) 

        noise_dict = {
            'l':noise[:,0],
            'l_pp': l_pp,
            'tt': tt_noise+tt_cv,
            'ee': ee_noise+ee_cv,
            'pp': pp_noise+pp_cv,
            'te': te_noise
        }

        fig,ax = plt.subplots()
        ax.plot(noise[:,0],tt_noise/tt_cv)
        ax.set_xlabel(r'$\ell$')
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_ylabel(r'$N_\ell^{TT}/C_\ell^{CV}$')
        fig.savefig('TT_noise.pdf')
        fig,ax = plt.subplots()
        ax.plot(noise[:,0],ee_noise/ee_cv)
        ax.set_xlabel(r'$\ell$')
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_ylabel(r'$N_\ell^{EE}/C_\ell^{CV}$')
        fig.savefig('EE_noise.pdf')
        fig,ax = plt.subplots()
        ax.plot(l_pp,pp_noise/pp_cv)
        ax.set_xlabel(r'$\ell$')
        ax.set_yscale('log')
        ax.grid(True)
        ax.set_ylabel(r'$N_\ell^{\phi\phi}/C_\ell^{CV}$')
        fig.savefig('pp_noise.pdf')

        return(noise_dict)

    def line_plots(self,FULL_cls,NN_cls,spectrum,save_path, N=100,cosmic_variance=False):
        """
        This plot function is used to plot spectra using both CLs from ClassFULL and CLassNET.   
        """

        if N>len(FULL_cls['kk']):
            print("More samples selected than caluclated. Reduce the number of plotted samples to {}".format(len(FULL_cls['kk'])))
            N = len(FULL_cls['kk'])

        if cosmic_variance and (spectrum not in HAS_COSMIC_VARIANCE):
            return

        #cut data
        x,y_FULL,y_NN,y_factor = self.cut_datasets(FULL_cls,NN_cls,spectrum,N)
        fig,ax = plt.subplots(figsize=self.figsize)

        if cosmic_variance:
            if spectrum == 'te':
                # get tt and ee spectra for cosmic variance te calculation
                _,y_FULL_tt,y_NN_tt,y_factor_tt = self.cut_datasets(FULL_cls,NN_cls,'tt',1)
                _,y_FULL_ee,y_NN_ee,y_factor_ee = self.cut_datasets(FULL_cls,NN_cls,'ee',1)

                # calc cosmic variance
                cosmic_var_te = np.sqrt(1/(2*x[0]+1)) * np.sqrt( y_FULL[0]**2 + y_FULL_tt[0] * y_FULL_ee[0] ) 
                ax.fill_between(x[0],-y_factor*cosmic_var_te,y_factor*cosmic_var_te,color="C0", label="Cosmic Variance",alpha=0.3)
            else:
                ax.set_yscale('log')
                cv = self.calculate_cosmic_variance(x[0],spectrum)
                ax.fill_between(x[0],LOG_YLIMITS[spectrum],y_factor*cv*y_FULL[0],color="C0", label="Cosmic Variance",alpha=0.3)
                ax.set_ylim((LOG_YLIMITS[spectrum],np.max(y_FULL[0]*y_factor)*10))



        if spectrum in ['pk','pk_cb']:
            for i in range(N):
                x_instance = x[i][x[i]!=0]
                y_FULL_instance= y_FULL[i][x[i]!=0]
                y_NN_instance= y_NN[i][x[i]!=0]
                y_factor_instance = y_factor[x[i]!=0]
                if i==0:
                    if spectrum == 'pk':
                        _ = "\mathrm{m}"
                    else:
                        _ = "\mathrm{cb}"
                    ax.plot(x_instance,y_FULL_instance*y_factor_instance,c='green',label=r'$P_{'+_+',\mathrm{Full}}$',linewidth = self.linewidth)
                    ax.plot(x_instance,y_NN_instance*y_factor_instance,c='blue',label= r"$P_{"+_+",\mathrm{Net}}$",linewidth = self.linewidth)
                    ax.plot(x_instance,abs(y_FULL_instance-y_NN_instance)*y_factor_instance,c='red',label=r"$|P_{"+_+",\mathrm{Dif}}|$",linewidth = self.linewidth)
                else:
                    ax.plot(x_instance,y_FULL_instance*y_factor_instance,c='green',linewidth = self.linewidth)
                    ax.plot(x_instance,y_NN_instance*y_factor_instance,c='blue',linewidth = self.linewidth)
                    ax.plot(x_instance,abs(y_FULL_instance-y_NN_instance)*y_factor_instance,c='red',linewidth = self.linewidth)

        else:
            for i in range(N):
                if i==0:
                    ax.plot(x[i],y_FULL[i]*y_factor,c='green',label=r"$C_{\ell,\mathrm{Full}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                    ax.plot(x[i],y_NN[i]*y_factor,c='blue',label=r"$C_{\ell,\mathrm{Net}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                    if spectrum == 'te':
                        ax.plot(x[i],(y_FULL[i]-y_NN[i])*y_factor,c='red',label=r"$C_{\ell,\mathrm{Diff}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                    else:
                        ax.plot(x[i],abs(y_FULL[i]-y_NN[i])*y_factor,c='red',label=r"$|C_{\ell,\mathrm{Diff}}^{%(x)s}|$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                else:
                    ax.plot(x[i],y_FULL[i]*y_factor,c='green',linewidth = self.linewidth)
                    ax.plot(x[i],y_NN[i]*y_factor,c='blue',linewidth = self.linewidth)
                    ax.plot(x[i],abs(y_FULL[i]-y_NN[i])*y_factor,c='red',linewidth = self.linewidth)
                    if spectrum == 'te':
                        ax.plot(x[i],(y_FULL[i]-y_NN[i])*y_factor,c='red',linewidth = self.linewidth)
                    else:
                        ax.plot(x[i],abs(y_FULL[i]-y_NN[i])*y_factor,c='red',linewidth = self.linewidth)

        if spectrum in LOG_XAXIS:
            ax.set_xscale('log')
        ax.set_xlim([np.min(x[x!=0]),np.max(x)])

        if spectrum in LOG_AXIS:
            ax.set_yscale('log')

        if spectrum in SYMLOG_AXIS:
            ax.set_yscale('symlog', linthresh=SYMLOG_LIN_THRESHOLD[spectrum])


        if spectrum in ['pk','pk_cb']:
            ax.set_xlabel(r'$k$ [$h$/Mpc$^{-1}$]')
            if spectrum == 'pk':
                ax.set_ylabel(r"${%(x)s}(k,z)$ [$h^{-1}$/Mpc$^{3}$]"%{"x":LATEX_DICT[spectrum]})
            else:
                ax.set_ylabel(r"${%(x)s}(k,z)$ [$h^{-1}$/Mpc$^{3}$]"%{"x":LATEX_DICT[spectrum]})
        else:
            ax.set_xlabel(r'$\ell$')
            if spectrum in HAS_ELL2:
                ax.set_ylabel(r"$\ell^2(\ell + 1)^2C_{\ell}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]})
            elif spectrum in HAS_ELL:
                ax.set_ylabel(r"$\ell(\ell + 1)C_{\ell}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]})
            else:
                ax.set_ylabel(r"$C_{\ell}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]})


        ax.grid(True)
        if cosmic_variance:
            ax.legend(loc='upper left')
        else:
            ax.legend()


        fig.savefig(save_path, bbox_inches='tight')
        plt.close()

    def line_plots_difference(self,FULL_cls,NN_cls,spectrum,save_path, N=100, measure='abs', noise=False, chisq = None, parameter = None, para_list= None):
        """
        This plot function is used to plot the difference using both CLs from ClassFULL and CLassNET.    
        """

        if (spectrum in SKIP_REL_DIF) and (measure=='rel'):
            return

        if (spectrum not in HAS_NOISE) and (noise==True):
            return
        elif (spectrum in HAS_NOISE) and (noise==True):
            noise_dict = self.load_planck_noise(FULL_cls)

        if N>len(FULL_cls['kk']):
            print("More samples selected than caluclated. Reduce the number of plotted samples to {}".format(len(FULL_cls['kk'])))
            N = len(FULL_cls['kk'])

        #cut data
        x,y_FULL,y_NN,y_factor = self.cut_datasets(FULL_cls,NN_cls,spectrum,N)

        fig,ax = plt.subplots(figsize=self.figsize)
        for i in range(N):
            # first: calcuate difference
            if noise==True:
                if spectrum=='pp':
                    x=x[:,:1998] #Cut down to 1998 datapoints
                    dif = np.nan_to_num((y_FULL[i,:1998]-y_NN[i,:1998])/y_FULL[i,:1998]/noise_dict[spectrum])
                elif spectrum =='te':
                    #We need tt and ee spectra for this
                    _,tt_FULL,_,_ = self.cut_datasets(FULL_cls,NN_cls,'tt',N)
                    _,ee_FULL,_,_ = self.cut_datasets(FULL_cls,NN_cls,'ee',N)

                    # calculate TE noise + variance
                    noise_te = np.sqrt(2/(2*x[i]+1)) * np.sqrt( y_FULL[i] * y_FULL[i] + noise_dict['tt'] * tt_FULL[i] * noise_dict['ee'] * ee_FULL[i] )

                    dif = np.nan_to_num((y_FULL[i]-y_NN[i])/noise_te)
                else:
                    dif = np.nan_to_num((y_FULL[i]-y_NN[i])/y_FULL[i]/noise_dict[spectrum])
            else:
                if spectrum in ['pk','pk_cb']:
                    x_instance = x[i][x[i]!=0]
                    y_FULL_instance= y_FULL[i][x[i]!=0]
                    y_NN_instance= y_NN[i][x[i]!=0]
                    y_factor_instance = y_factor[x[i]!=0]
                    if measure == 'abs':
                        dif = (y_FULL_instance-y_NN_instance)*y_factor_instance
                    elif measure == 'rel':
                        dif = np.nan_to_num((y_FULL_instance-y_NN_instance)/y_FULL_instance)
                else:
                    if measure == 'abs':
                        dif = (y_FULL[i]-y_NN[i])*y_factor
                    elif measure == 'rel':
                        dif = np.nan_to_num((y_FULL[i]-y_NN[i])/y_FULL[i])
                    else:
                        raise ValueError('measure either abs or rel')
            
            # second: plot 
            if chisq is None:
                if parameter is not None:
                    # cut parameter list if necessary
                    para_list_cut = para_list[:max(N,len(para_list))]

                    # normalize parameter values
                    max_val = np.max(para_list_cut)
                    min_val = np.min(para_list_cut)

                    norm_val = (para_list_cut[i] - min_val) / (max_val - min_val)

                    color = plt.cm.plasma(norm_val)
                    if spectrum in ['pk','pk_cb']:
                        ax.plot(x_instance,dif,c=color,linewidth = self.linewidth)
                    else:
                        ax.plot(x[i],dif,c=color,linewidth = self.linewidth)

                else:
                    if spectrum in ['pk','pk_cb']:
                        ax.plot(x_instance,dif,c='red',linewidth = self.linewidth)
                    else:
                        ax.plot(x[i],dif,c='red',linewidth = self.linewidth)
            else:
                # constants for the colorbar
                max_sigma = 5
                scale_power = 3

                # Here we calculate from chisquare to the distance to the center of the ellipsoid in sigmata.
                prob = sc.stats.chi2.sf(chisq[i],9) #here we make an asumption on the dof. To be changed in other cases.
                n_sigma = -sc.stats.norm.ppf(prob/2)

                norm_sigma = n_sigma**scale_power / max_sigma**scale_power

                color = plt.cm.gist_rainbow(norm_sigma)
                if spectrum in ['pk','pk_cb']:
                    ax.plot(x_instance,dif,c=color,linewidth = self.linewidth)
                else:
                    ax.plot(x[i],dif,c=color,linewidth = self.linewidth)

        if spectrum in LOG_XAXIS:
            ax.set_xscale('log')
        ax.set_xlim([np.min(x[x!=0]),np.max(x)])

        if chisq is not None:
            cmap = mpl.cm.gist_rainbow
            norm = mpl.colors.PowerNorm(gamma=scale_power, vmin=0.0,vmax=max_sigma)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm,boundaries=np.arange(-0.05,2.1,.1), label=r'distance from center of domain [$\sigma$]', ticks=np.arange(5)+1)         
        
        if parameter is not None:
            cmap = mpl.cm.plasma
            norm = mpl.colors.Normalize(vmin=min_val,vmax=max_val)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            import matplotlib.ticker as ticker
            cb = plt.colorbar(sm,boundaries=np.arange(-0.05,2.1,.1), label=r"$%(x)s$"%{"x":LATEX_DICT_PARAMETER[parameter]})   
            cb.locator = ticker.MaxNLocator(10)
            cb.update_ticks()


        if measure == 'abs':
            if spectrum in ['pk','pk_cb']:
                if spectrum == 'pk':
                    _ = '\mathrm{m}'
                else:
                    _ = '\mathrm{cb}'
                ax.set_xlabel(r'$k$ [$h$/Mpc$^{-1}$]')
                ax.set_ylabel(r'$(P_{'+_+',\mathrm{Full}}(k,z)-P_{'+_+',\mathrm{Net}}(k,z))$ [$h^{-1}$/Mpc$^{3}$]')
            else:
                ax.set_xlabel(r'$\ell$')
                if spectrum in HAS_ELL2:
                    ax.set_ylabel(r"$\ell^2(\ell + 1)^2(C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s})$"%{"x":LATEX_DICT[spectrum]})
                elif spectrum in HAS_ELL:
                    ax.set_ylabel(r"$\ell(\ell + 1)(C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s})$"%{"x":LATEX_DICT[spectrum]})
                else:
                    ax.set_ylabel(r"$(C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s})$"%{"x":LATEX_DICT[spectrum]})
            # if spectrum in SYMLOG_AXIS:
            #     plt.yscale('symlog', linthresh=SYMLOG_LIN_THRESHOLD[spectrum])
        elif measure == 'rel':
            if spectrum in ['pk','pk_cb']:
                if spectrum == 'pk':
                    _ = '\mathrm{m}'
                else:
                    _ = '\mathrm{cb}'
                ax.set_xlabel(r'$k$ [$h$/Mpc$^{-1}$]')
                ax.set_ylabel(r'$(P_{'+_+',\mathrm{Full}}(k,z)-P_{'+_+',\mathrm{Net}}(k,z))/P_{'+_+',\mathrm{Full}}(k,z)$ ')
            else:
                ax.set_xlabel(r'$\ell$')
                ax.set_ylabel(r"$(C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s})/C_{\ell,\mathrm{Full}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]})

        #if spectrum in LOG_AXIS:
        #    ax.set_yscale('log')

        ax.grid(True)

        fig.savefig(save_path, bbox_inches='tight')
        plt.close()

    def line_plots_contain(self,FULL_cls,NN_cls,spectrum,save_path, N=1000, sigmas=[0.68,0.95,0.99], noise = False):
        """
        This plot function is used to plot the difference using both CLs from ClassFULL and CLassNET.    
        """
        if N>len(FULL_cls['kk']):
            print("More samples selected than caluclated. Reduce the number of plotted samples to {}".format(len(FULL_cls['kk'])))
            N = len(FULL_cls['kk'])
        
        # Check for noise
        if (spectrum not in HAS_NOISE) and (noise==True):
            return
        elif (spectrum in HAS_NOISE) and (noise==True):
            noise_dict = self.load_planck_noise(FULL_cls)

        #cut data
        x,y_FULL,y_NN,y_factor = self.cut_datasets(FULL_cls,NN_cls,spectrum,N)
        
        #define difference
        if spectrum in ['pk','pk_cb']:
            y_difference = np.zeros((N,700))
        else:
            y_difference = np.zeros((N,len(x[0])))

        if noise:
            if spectrum not in HAS_NOISE:
                return 
            else:
                if spectrum=='pp':
                    x=x[:,:1998] #Cut down to 1998 datapoints
                    y_difference = np.zeros((N,len(x[0])))
                    for i in range(N):
                        y_difference[i] = np.nan_to_num(abs(y_FULL[i,:1998]-y_NN[i,:1998])/y_FULL[i,:1998]/noise_dict[spectrum])

                elif spectrum =='te':
                    for i in range(N):
                        y_difference[i] = np.nan_to_num(abs(y_FULL[i]-y_NN[i])/noise_dict[spectrum])

                else:
                    for i in range(N):
                        y_difference[i] = np.nan_to_num(abs(y_FULL[i]-y_NN[i])/y_FULL[i]/noise_dict[spectrum])
        else:
            if spectrum in SKIP_REL_DIF:
                for i in range(N):
                    y_difference[i] = abs(y_FULL[i]-y_NN[i])*y_factor
            elif spectrum in ['pk','pk_cb']:
                k_min = x[:,0]
                x_new = np.repeat([np.logspace(-4,1,num=700)],2,axis=0)

                for i in range(N):
                    # for pk we need to interpolate the values to get the same k bins
                    full = sc.interpolate.interp1d(x[i], y_FULL[i])
                    y_FULL_interpolate = full(x_new[0])
                    net = sc.interpolate.interp1d(x[i], y_NN[i])
                    y_NN_interpolate = net(x_new[0])
                    y_difference[i] = np.nan_to_num(abs(y_FULL_interpolate-y_NN_interpolate)/y_FULL_interpolate)
                x = x_new

            else:
                for i in range(N):
                    y_difference[i] = np.nan_to_num(abs(y_FULL[i]-y_NN[i])/y_FULL[i])

        fig,ax = plt.subplots(figsize=self.figsize)

        # calculate sigma curves here
        sigma_curves = np.zeros((len(x[0]),len(sigmas)))

        for i,val in enumerate(x[0]):
            if spectrum in ['pk','pk_cb']:
                for j,sigma in enumerate(sigmas):
                    mask = (k_min<val)
                    sigma_curves[i,j] = np.sort(y_difference[mask,i])[int(sum(mask)*sigma)-1]
            else:
                for j,sigma in enumerate(sigmas):
                    sigma_curves[i,j] = np.sort(y_difference[:,i])[int(N*sigma)-1]

        for j,sigma in enumerate(sigmas):
            ax.fill_between(x[0],sigma_curves[:,j],color='red', label = str(int(sigma*100)) + ' %', alpha=0.2+0.2*j)
        
        if (spectrum in LOG_AXIS) and (spectrum not in ['pk','pp']):
            ax.set_yscale('log')
            ax.set_ylim([np.min(sigma_curves)*0.5,np.max(sigma_curves)*1.2])
        else:
            ax.set_ylim([0,np.max(sigma_curves)*1.2])

        if spectrum in LOG_XAXIS:
            ax.set_xscale('log')

        ax.set_xlim([np.min(x[x!=0]),np.max(x)])


        if noise:
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r"$|C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s}|/\sigma_{\ell}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]})
        else:
            if spectrum in SKIP_REL_DIF:
                ax.set_xlabel(r'$\ell$')
                ax.set_ylabel(r"$\ell(\ell + 1)|C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s}|$"%{"x":LATEX_DICT[spectrum]})
            else:
                if spectrum in ['pk','pk_cb']:
                    ax.set_ylabel(r'$|P_{\mathrm{Full}}(k,z)-P_{\mathrm{Net}}(k,z)|/P_{\mathrm{Full}}(k,z)$')
                    ax.set_xlabel(r'$k$ [$h$/Mpc$^{-1}$]')
                    if spectrum == 'pk':
                        _name_full = 'P_{\mathrm{m,Full}}'
                        _name_net = 'P_{\mathrm{m,Full}}'
                    else:
                        _name_full = 'P_{\mathrm{cb,Full}}'
                        _name_net = 'P_{\mathrm{cb,Full}}'

                    ax.set_ylabel(r"$|%(x)s-%(y)s|/%(x)s$"%{"x":_name_full, "y":_name_net})


                else:
                    ax.set_xlabel(r'$\ell$')
                    ax.set_ylabel(r"$|C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s}|/|C_{\ell,\mathrm{Full}}^{%(x)s}|$"%{"x":LATEX_DICT[spectrum]})


        ax.legend()
        ax.grid(True)

        fig.savefig(save_path, bbox_inches='tight')
        plt.close()


