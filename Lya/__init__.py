from montepython.likelihood_class import Likelihood
import os
import numpy as np
#from scipy.interpolate import InterpolatedUnivariateSpline
from lmfit import Minimizer, Parameters, report_fit
import matplotlib.pyplot as plt
import time

class Lya(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        self.need_cosmo_arguments(data, {'output': 'mPk'})
        self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': 1.5*self.khmax}) 
        #self.need_cosmo_arguments(data, {'z_max_pk': zmax})

        #self.z = np.linspace(self.zmin, self.zmax, num=self.nbin)
        self.kh_0 = np.linspace(self.khmin, self.khmax, num=self.k_size)

        self.file_exists = False
        file_path = os.path.join(self.data_directory, self.grid_file)
        if os.path.exists(file_path):
           self.file_exists = True

#        if self.file_exist is True:
#           with open(file_path, 'r') as grid_file:
#                line = grid_file.readline()
#                while line.find('#') != -1:
#                    line = grid_file.readline()
#                while (line.find('\n') != -1 and len(line) == 1):
#                    line = grid_file.readline()
#                for index in xrange(self.grid_size):
#                    self.tau_grid[index] = float(line.split()[0])#check column
#                    self.sigma8_grid[index] = float(line.split()[0])
#                    self.neff_grid[index] = float(line.split()[0])
#                    self.alpha[index] = float(line.split()[0])
#                    self.beta[index] = float(line.split()[0])
#                    self.gamma[index] = float(line.split()[0])
#                    line = grid_file.readline()
#        else
#           raise io_mp.ConfigurationError('Error: grid file is missing')
#           exit()
        return

    def loglkl(self, cosmo, data):

        #see likelihood_class get_flat_fid
        print data.cosmo_arguments
        param_backup = data.cosmo_arguments
        #pba->Omega0_idr = pba->stat_f_idr*pow(pba->xi_idr,4.)*pba->Omega0_g;
        #N_dark = pba->Omega0_idr/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;
        DeltaNeff=data.cosmo_arguments['stat_f_idr']*(data.cosmo_arguments['csi_idr']**4)/7.*8./((4./11.)**(4./3.))
        eta2=(1.+0.2271*(data.cosmo_arguments['N_ur']+DeltaNeff))/(1.+0.2271*data.cosmo_arguments['N_ur'])
        eta=np.sqrt(eta2)
        ob=eta2*data.cosmo_arguments['omega_b']
        oc=eta2*data.cosmo_arguments['omega_cdm']
        h=eta*data.cosmo_arguments['h']
        print DeltaNeff,eta2,ob,oc,h
        #convert Neff according to arXiv:
        data.cosmo_arguments = {'output': ' mPk ','P_k_max_h/Mpc': 1.5*self.khmax,
                                'omega_b': ob,'omega_cdm': oc,'h':h,
                                'A_s': self.cosmo_arguments['A_s'],'n_s': self.cosmo_arguments['n_s'], 'tau_reio': self.cosmo_arguments['tau_reio'],
                                'N_ur':'3.046'}
        print data.cosmo_arguments
        cosmo.empty()
        cosmo.set(data.cosmo_arguments)
        cosmo.compute(['lensing'])

        Plin_equiv_0 = np.zeros(len(kh_0), 'float64')
        for index_k in range(len(self.kh_0)):
            Plin_equiv_0[index_k] = cosmo.pk_lin(self.kh_0[index_k]*h, 0.0)
        Plin_equiv_0 *= h**3
        z_reio=cosmo.z_reio()
        neff=cosmo.neff()
        sigma8=cosmo.sigma8()
        print z_reio, sigma8, neff

        cosmo.empty()
        data.cosmo_arguments = param_backup
        h=data.cosmo_arguments['h']
        print data.cosmo_arguments
        cosmo.set(data.cosmo_arguments)
        cosmo.compute(['lensing'])

        Plin_0 = np.zeros(len(kh_0), 'float64')
        for index_k in range(len(self.kh_0)):
            Plin_0[index_k] = cosmo.pk_lin(self.kh_0[index_k]*h, 0.0)
        Plin_0 *= h**3

        Tk_0 = np.sqrt(Plin_0/Plin_equiv_0) 

        #Now merge with Riccardo's interpolation code
        #setting k_max (i.e. cutting oscillations from the fitted region)
        for index_k in range(len(kh_0)):
            index_khmax = -1
            if Plin_0[index_k] <= 0.01: #and any(data[i:]) > 1.05*data[i]: #MArchi perhaps this criteria has to be adjusted?
               index_khmax = index_k
               print index_khmax
               break

        kh = kh_0[:index_khmax]
        Plin_equiv = Plin_equiv_0[:index_khmax]
        Plin = Plin_0[:index_khmax]

        Tk = np.sqrt(Plin/Plin_equiv)

        # fitting the given linear P(k) with the {alpha,beta,gamma}-formula 

        #model function #T^2=P_model/P_ref
        def T(kh,alpha,beta,gamma):
            return (1. + (alpha*kh)**(beta))**(gamma)

        # define objective function: returns the array to be minimized
        def fcn2min(params, kh, Tk):
            alpha = params['alpha']
            beta = params['beta']
            gamma = params['gamma']
            model = T(kh,alpha,beta,gamma) #(1. + (alpha*kappa_interp)**(beta))**(gamma)
            return (model - Tk)      #standard residuals

        # create a set of Parameters
        params = Parameters()
        params.add('alpha', value=0.001, min = 0., max = 0.3)
        params.add('beta', value=2.24, min = 0.5, max = 10.)
        params.add('gamma', value=-4.46, min=-10., max=-0.1)

        # do fit, default is with least squares method
        t0_fit = time.clock()

        minner = Minimizer(fcn2min, params, fcn_args=(kh, Tk))
        result = minner.minimize(method = 'leastsq')
        best_alpha = result.params['alpha'].value
        best_beta  = result.params['beta'].value
        best_gamma = result.params['gamma'].value

        t1_fit = time.clock()

        # write error report
        report_fit(result)

        plt.xlabel('k [h/Mpc]')
        plt.ylabel('$P_{nCDM}/P_{CDM}$')

        plt.ylim(0.,1.1)
        plt.xlim(0.1,150.)
        plt.xscale('log')
        #plt.yscale('log')
        plt.grid(True)

        plt.plot(kh_0, Tk_0**2, 'r')
        plt.plot(kh_0, (T(k_0, best_alpha, best_beta, best_gamma))**2, 'b--')
        #plt.show()
        plt.savefig('grid_fit_plot.pdf')


        return
