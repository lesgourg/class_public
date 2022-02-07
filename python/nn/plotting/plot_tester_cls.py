import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import pkg_resources

# Settings for different spectra
REQUIERES_ELL_FACTOR = ['tt','bb','ee','te','tp']
SKIP_REL_DIF = ['te','tp']
LOG_AXIS = ['pk','pp','tp']
LOG_XAXIS = ['ee','te','pk']
HAS_COSMIC_VARIANCE = ['tt','pp','ee']
HAS_NOISE = ['pp','tt','ee','te']
LOG_YLIMITS= {'tt':1e-13,'pp':1e-18,'ee':1e-16}

# latex dict
LATEX_DICT = {
    'tt':r'\mathrm{TT}',
    'bb':r'\mathrm{BB}',
    'ee':r'\mathrm{EE}',
    'te':r'\mathrm{TE}',
    'pp':r'\phi \phi',
    'pk':r'P_k',
    'tp':r'\mathrm{T} \phi'
}
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
        if spectrum == 'pk':
            # Load data
            x = np.array(FULL_cls['kk'])
            y_FULL = np.array(FULL_cls['pk'])
            y_NN = np.array(NN_cls['pk'])

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

            y_factor = x[0] * (x[0] + 1)

        return(x,y_FULL,y_NN,y_factor)

    def calculate_cosmic_variance(self, ell, spectrum):

        if spectrum in ["ee","tt","pp"]:
            cv = np.sqrt(2/(2*ell+1))
            return cv

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
        pp_noise = np.sqrt(2/(2*l_pp+1))*noise[:1998,3]/FULL_cls['pp'][0][2:2000]/(2.726*10**6)**2

        # caclualte TE component
        cv_term = FULL_cls['te'][0][2:] * FULL_cls['te'][0][2:]
        noise_cv_term = (FULL_cls["tt"][0][2:] + np.array(noise[:,1])/(2.726*10**6)**2)  *  (FULL_cls["ee"][0][2:] + np.array(noise[:,2])/(2.726*10**6)**2)
        te_noise = np.sqrt(1/(2*FULL_cls["ell"][0][2:]+1)) * np.sqrt(cv_term + noise_cv_term) 

        noise_dict = {
            'l':noise[:,0],
            'l_pp': l_pp,
            'tt': tt_noise+tt_cv,
            'ee': ee_noise+ee_cv,
            'pp': pp_noise*l_pp*(l_pp+1)/(2*np.math.pi)*l_pp*(l_pp+1)+pp_cv,
            'te': te_noise
        }
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
            ax.set_yscale('log')
            cv = self.calculate_cosmic_variance(x[0],spectrum)
            ax.fill_between(x[0],LOG_YLIMITS[spectrum],y_factor*cv*y_FULL[0],color="C0", label="Cosmic Variance",alpha=0.3)
            ax.set_ylim((LOG_YLIMITS[spectrum],np.max(y_FULL[0]*y_factor)*10))

        if spectrum=='pk':
            for i in range(N):
                x_instance = x[i][x[i]!=0]
                y_FULL_instance= y_FULL[i][x[i]!=0]
                y_NN_instance= y_NN[i][x[i]!=0]
                y_factor_instance = y_factor[x[i]!=0]
                if i==0:
                    ax.plot(x_instance,y_FULL_instance*y_factor_instance,c='green',label=r"$P_{k,\mathrm{Full}}$",linewidth = self.linewidth)
                    ax.plot(x_instance,y_NN_instance*y_factor_instance,c='blue',label=r"$P_{k,\mathrm{Net}}$",linewidth = self.linewidth)
                    ax.plot(x_instance,abs(y_FULL_instance-y_NN_instance)*y_factor_instance,c='red',label=r"$|P_{k,\mathrm{Dif}}|$",linewidth = self.linewidth)
                else:
                    ax.plot(x_instance,y_FULL_instance*y_factor_instance,c='green',linewidth = self.linewidth)
                    ax.plot(x_instance,y_NN_instance*y_factor_instance,c='blue',linewidth = self.linewidth)
                    ax.plot(x_instance,abs(y_FULL_instance-y_NN_instance)*y_factor_instance,c='red',linewidth = self.linewidth)

        else:
            for i in range(N):
                if i==0:
                    ax.plot(x[i],y_FULL[i]*y_factor,c='green',label=r"$C_{\ell,\mathrm{Full}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                    ax.plot(x[i],y_NN[i]*y_factor,c='blue',label=r"$C_{\ell,\mathrm{Net}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                    ax.plot(x[i],abs(y_FULL[i]-y_NN[i])*y_factor,c='red',label=r"$|C_{\ell,\mathrm{Diff}}^{%(x)s}|$"%{"x":LATEX_DICT[spectrum]},linewidth = self.linewidth)
                else:
                    ax.plot(x[i],y_FULL[i]*y_factor,c='green',linewidth = self.linewidth)
                    ax.plot(x[i],y_NN[i]*y_factor,c='blue',linewidth = self.linewidth)
                    ax.plot(x[i],abs(y_FULL[i]-y_NN[i])*y_factor,c='red',linewidth = self.linewidth)

        if spectrum in LOG_XAXIS:
            ax.set_xscale('log')
        ax.set_xlim([np.min(x[x!=0]),np.max(x)])

        if spectrum in LOG_AXIS:
            ax.set_yscale('log')


        if spectrum=='pk':
            ax.set_xlabel(r'$k$')
            ax.set_ylabel(r'$P_k$')
        else:
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r"$\ell(\ell + 1)(C_{\ell}^{%(x)s})$"%{"x":LATEX_DICT[spectrum]})


        ax.grid(True)
        ax.legend()

        fig.savefig(save_path, bbox_inches='tight')
        del fig




    def line_plots_difference(self,FULL_cls,NN_cls,spectrum,save_path, N=100, measure='abs', noise=False, chisq = None):
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
                if spectrum=='pk':
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
                        if spectrum=='tt':
                            if dif[100]>0.1e-10:
                                print(i)
                    elif measure == 'rel':
                        dif = np.nan_to_num((y_FULL[i]-y_NN[i])/y_FULL[i])
                    else:
                        raise ValueError('measure either abs or rel')
            
            if chisq is None:
                if spectrum=='pk':
                    ax.plot(x_instance,dif,c='red',linewidth = self.linewidth)
                else:
                    ax.plot(x[i],dif,c='red',linewidth = self.linewidth)
            else:
                max_chi = np.max(chisq)
                norm_chisq = chisq[i] / max_chi
                color = plt.cm.jet(norm_chisq)
                if spectrum=='pk':
                    ax.plot(x_instance,dif,c=color,linewidth = self.linewidth)
                else:
                    ax.plot(x[i],dif,c=color,linewidth = self.linewidth)

        if spectrum in LOG_XAXIS:
            ax.set_xscale('log')
        ax.set_xlim([np.min(x[x!=0]),np.max(x)])

        if chisq is not None:
            cmap = mpl.cm.jet
            norm = mpl.colors.Normalize(vmin=0,vmax=max_chi)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm,boundaries=np.arange(-0.05,2.1,.1), label='chisq', ticks=4*np.arange(int(max_chi/4)+1))         

        if measure == 'abs':
            if spectrum=='pk':
                ax.set_xlabel(r'$k$')
                ax.set_ylabel(r'$P_{k,\mathrm{Full}}-P_{k,\mathrm{Net}}$')
            else:
                ax.set_xlabel(r'$\ell$')
                ax.set_ylabel(r"$\ell(\ell + 1)(C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s})$"%{"x":LATEX_DICT[spectrum]})
        elif measure == 'rel':
            if spectrum=='pk':
                ax.set_xlabel(r'$k$')
                ax.set_ylabel(r'$(P_{k,\mathrm{Full}}-P_{k,\mathrm{Net}})/P_{k,\mathrm{Full}}$')
            else:
                ax.set_xlabel(r'$\ell$')
                ax.set_ylabel(r"$(C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s})/C_{\ell,\mathrm{Full}}^{%(x)s}$"%{"x":LATEX_DICT[spectrum]})

        #if spectrum in LOG_AXIS:
        #    ax.set_yscale('log')

        ax.grid(True)

        fig.savefig(save_path, bbox_inches='tight')

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

        # Containment plots for pk is difficult due to different coverage of k values and will be skipped here. Maybe something for somewhen else
        if spectrum=='pk':
            return

        #cut data
        x,y_FULL,y_NN,y_factor = self.cut_datasets(FULL_cls,NN_cls,spectrum,N)
        
        #define difference
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
            else:
                for i in range(N):
                    y_difference[i] = np.nan_to_num(abs(y_FULL[i]-y_NN[i])/y_FULL[i])

        fig,ax = plt.subplots(figsize=self.figsize)

        # calculate sigma curves here
        sigma_curves = np.zeros((len(x[0]),len(sigmas)))

        for i in range(len(x[0])):
            for j,sigma in enumerate(sigmas):
                sigma_curves[i,j] = np.sort(y_difference[:,i])[int(N*sigma)-1]

        for j,sigma in enumerate(sigmas):
            ax.fill_between(x[0],sigma_curves[:,j],color='red', label = str(int(sigma*100)) + ' %', alpha=0.2+0.2*j)
        
        if spectrum in LOG_AXIS:
            ax.set_yscale('log')
            ax.set_ylim([np.min(sigma_curves),np.max(sigma_curves)*1.2])
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
                if spectrum=='pk':
                    ax.set_xlabel(r'$k$')
                    ax.set_ylabel(r'$|P_{k,\mathrm{Full}}-P_{k,\mathrm{Net}}|$')
                else:
                    ax.set_xlabel(r'$\ell$')
                    ax.set_ylabel(r"$\ell(\ell + 1)|C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s}|$"%{"x":LATEX_DICT[spectrum]})
            else:
                ax.set_xlabel(r'$\ell$')
                ax.set_ylabel(r"$|C_{\ell,\mathrm{Full}}^{%(x)s}-C_{\ell,\mathrm{Net}}^{%(x)s}|/|C_{\ell,\mathrm{Full}}^{%(x)s}|$"%{"x":LATEX_DICT[spectrum]})


        ax.legend()
        ax.grid(True)

        fig.savefig(save_path, bbox_inches='tight')


