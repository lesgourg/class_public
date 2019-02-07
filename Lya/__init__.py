from montepython.likelihood_class import Likelihood
import io_mp
import re  # Module to handle regular expressions
#from datetime import date
import sys
import os
import numpy as np
#from copy import deepcopy
from scipy import interpolate
from lmfit import Minimizer, Parameters, report_fit
from scipy.linalg import block_diag
import pprint, pickle
#import matplotlib.pyplot as plt
#import time

#Lyman alpha likelihood by M. Archidiacono, R. Murgia, D.C. Hooper, J. Lesgourgues, M. Viel

class Lya(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        print "Initializing Lya likelihood"

        self.need_cosmo_arguments(data, {'output': 'mPk'})
        self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': 1.5*self.kmax})
        self.need_cosmo_arguments(data, {'compute_neff_Lya':'yes'})
        self.need_cosmo_arguments(data, {'Lya_k_s_over_km':self.k_s_over_km})
        self.need_cosmo_arguments(data, {'Lya_z':self.z})

        #lcdm_points = 33    #number of grid points for the lcdm case (i.e. alpha=0, regardless of beta and gamma values)
        self.params_numbers = 3  #number of non-astro params (i.e. alpha,beta and gamma)

        alphas = np.zeros(self.grid_size, 'float64')
        betas = np.zeros(self.grid_size, 'float64')
        gammas = np.zeros(self.grid_size, 'float64')

        self.bin_file_path = os.path.join(command_line.folder,self.bin_file_name)
        if not os.path.exists(self.bin_file_path):
           with open(self.bin_file_path, 'w') as bin_file:
                bin_file.write('#')
                for name in data.get_mcmc_parameters(['varying']):
                    name = re.sub('[$*&]', '', name)
                    bin_file.write(' %s\t' % name)
                for name in data.get_mcmc_parameters(['derived']):
                    name = re.sub('[$*&]', '', name)
                    bin_file.write(' %s\t' % name)
                bin_file.write('\n')
                bin_file.close()
        if 'z_reio' not in data.get_mcmc_parameters(['derived']) or 'sigma8' not in data.get_mcmc_parameters(['derived']) or 'neff' not in data.get_mcmc_parameters(['derived']):
           raise io_mp.ConfigurationError('Error: Lya likelihood need z_reio, sigma8 and neff as derived parameters')
           exit()

        file_path = os.path.join(self.data_directory, self.grid_file)
        if os.path.exists(file_path):
           with open(file_path, 'r') as grid_file:
                line = grid_file.readline()
                while line.find('#') != -1:
                    line = grid_file.readline()
                while (line.find('\n') != -1 and len(line) == 3):
                    line = grid_file.readline()
                for index in xrange(self.grid_size):
                    alphas[index] = float(line.split()[0])
                    betas[index] = float(line.split()[1])
                    gammas[index] = float(line.split()[2])
                    line = grid_file.readline()
                grid_file.close()
        else:
           raise io_mp.ConfigurationError('Error: grid file is missing')
           exit()

        X_real = np.zeros((self.grid_size, self.params_numbers),'float64') #real params

        for k in range(self.grid_size):  #real params
           X_real[k][0] = self.khalf(alphas[k], betas[k], gammas[k])    #k_1/2
           X_real[k][1] = betas[k]
           X_real[k][2] = gammas[k]

        #for the normalization (see Alex's notes) #???
        self.a_min = min(X_real[:,0])
        self.b_min = min(X_real[:,1])
        self.g_min = min(X_real[:,2])
        self.a_max = max(X_real[:,0])
        self.b_max = max(X_real[:,1])
        self.g_max = max(X_real[:,2])

        #redshift independent parameters - params order: z_reio, sigma_8, n_eff, f_UV
        self.zind_param_size = [3, 5, 5, 3] #how many values I have for each param
        self.zind_param_min = np.array([7., 0.5, -2.6, 0.])
        self.zind_param_max = np.array([15., 1.5, -2.0, 1.])
        zind_param_ref = np.array([9., 0.829, -2.3074, 0.])
        self.zreio_range = self.zind_param_max[0]-self.zind_param_min[0]
        self.neff_range = self.zind_param_max[2]-self.zind_param_min[2]

        # redshift dependent parameters - params order: params order: mean_f , t0, slope
        zdep_params_size = [9, 3, 3] #how many values I have for each param
        zdep_params_refpos = [4, 1, 2] #where to store the P_F(ref)# DATA

        #MEAN FLUXES values###
        flux_ref_old = (np.array([0.669181, 0.617042, 0.564612, 0.512514, 0.461362, 0.411733, 0.364155, 0.253828, 0.146033, 0.0712724]))
        #flux_min_meanf = (np.array([0.401509, 0.370225, 0.338767, 0.307509, 0.276817, 0.24704, 0.218493, 0.152297, 0.0876197, 0.0427634]))
        #flux_max_meanf = (np.array([0.936854, 0.863859, 0.790456, 0.71752, 0.645907, 0.576426, 0.509816, 0.355359, 0.204446, 0.0997813]))

        # DATA
        # FIRST (NOT USED) DATASET (19 wavenumbers) ***XQ-100***
        self.zeta_range_XQ = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2]  #list of redshifts corresponding to the 19 wavenumbers (k) (see ./SPECTRA/README)
        self.k_XQ = [0.003,0.006,0.009,0.012,0.015,0.018,0.021,0.024,0.027,0.03,0.033,0.036,0.039,0.042,0.045,0.048,0.051,0.054,0.057]

        # SECOND DATASET (7 wavenumbers) ***HIRES/MIKE***
        self.zeta_range_mh = [4.2, 4.6, 5.0, 5.4]  #list of redshifts corresponding to the 7 wavenumbers (k) (see ./SPECTRA/README)
        self.k_mh = [0.00501187,0.00794328,0.0125893,0.0199526,0.0316228,0.0501187,0.0794328] #in s/km

        self.zeta_full_length = (len(self.zeta_range_XQ) + len(self.zeta_range_mh))
        self.kappa_full_length = (len(self.k_XQ) + len(self.k_mh))
        self.redshift = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.2, 4.6, 5.0, 5.4] #which snapshots (first 7 for first dataset, last 4 for second one)

        #T0 AND SLOPE VALUES
        t0_ref_old = np.array([11251.5, 11293.6, 11229.0, 10944.6, 10421.8, 9934.49, 9227.31, 8270.68, 7890.68, 7959.4])
        slope_ref_old = np.array([1.53919, 1.52894, 1.51756, 1.50382, 1.48922, 1.47706, 1.46909, 1.48025, 1.50814, 1.52578])

        t0_values_old = np.zeros(( 10, zdep_params_size[1] ),'float64')
        t0_values_old[:,0] = np.array([7522.4, 7512.0, 7428.1, 7193.32, 6815.25, 6480.96, 6029.94, 5501.17, 5343.59, 5423.34])
        t0_values_old[:,1] = t0_ref_old[:]
        t0_values_old[:,2] = np.array([14990.1, 15089.6, 15063.4, 14759.3, 14136.3, 13526.2, 12581.2, 11164.9, 10479.4, 10462.6])

        slope_values_old = np.zeros(( 10, zdep_params_size[2] ),'float64')
        slope_values_old[:,0] = np.array([0.996715, 0.979594, 0.960804, 0.938975, 0.915208, 0.89345, 0.877893, 0.8884, 0.937664, 0.970259])
        slope_values_old[:,1] = [1.32706, 1.31447, 1.30014, 1.28335, 1.26545, 1.24965, 1.2392, 1.25092, 1.28657, 1.30854]
        slope_values_old[:,2] = slope_ref_old[:]

        self.t0_min = t0_values_old[:,0]*0.1
        self.t0_max = t0_values_old[:,2]*1.4
        self.slope_min = slope_values_old[:,0]*0.8
        self.slope_max = slope_values_old[:,2]*1.15

        #IMPORTING THE TWO GRIDS FOR KRIGING
        # Here I import the grids that I pre-computed through the file "setting_Kriging_grid_2R.py"
        file_path = os.path.join(self.data_directory, self.astro_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           self.input_full_matrix_interpolated_ASTRO = pickle.load(pkl)
           #print self.full_matrix_interpolated_ASTRO.shape
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: astro spectra file is missing')
           exit()

        file_path = os.path.join(self.data_directory, self.abg_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           self.input_full_matrix_interpolated_ABG = pickle.load(pkl)
           #print self.full_matrix_interpolated_ABG.shape
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: abg spectra file is missing')
           exit()

        ALL_zdep_params = len(flux_ref_old) + len(t0_ref_old) + len(slope_ref_old)
        grid_length_ABG = len(self.input_full_matrix_interpolated_ABG[0,0,:])
        grid_length_ASTRO = len(self.input_full_matrix_interpolated_ASTRO[0,0,:])
        astroparams_number_KRIG = len(self.zind_param_size) + ALL_zdep_params

        #### --- ABG GRID --- {alpha, beta, gamma} GRID
        file_path = os.path.join(self.data_directory, self.abg_grid_file)
        if os.path.exists(file_path):
           self.X_ABG = np.zeros((grid_length_ABG, self.params_numbers), 'float64')
           for param_index in range(self.params_numbers):
               self.X_ABG[:,param_index] = np.genfromtxt(file_path, usecols=[param_index], skip_header=1)
        else:
           raise io_mp.ConfigurationError('Error: abg grid file is missing')
           exit()

        #### --- ASTRO GRID --- ORDER OF THE COMPLETE LIST OF ASTRO PARAMS: z_reio, sigma_8, n_eff, f_UV, mean_f(z), t0(z), slope(z)
        #### HIGH REDSHIFT
        file_path = os.path.join(self.data_directory, self.abg_astro_grid_file)
        if os.path.exists(file_path):
           self.X = np.zeros((grid_length_ASTRO,astroparams_number_KRIG), 'float64')
           for param_index in range(astroparams_number_KRIG):
               self.X[:,param_index] = np.genfromtxt(file_path, usecols=[param_index], skip_header=1)
        else:
           raise io_mp.ConfigurationError('Error: abg+astro grid file is missing')
           exit()

        #  STUFF FOR INTERPOLATING IN THE ASTROPARAMS SPACE ####################################
        self.redshift_list = np.array([3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4]) #combined dataset (MIKE/HIRES + XQ-100)
        #astrokrig_result = np.zeros((self.zeta_full_length, self.kappa_full_length ), 'float64')
        #minimum and maximum values for the kriging normalisation###
        self.F_prior_min = np.array([0.535345,0.493634,0.44921,0.392273,0.338578,0.28871,0.218493,0.146675,0.0676442,0.0247793])
        self.F_prior_max = np.array([0.803017,0.748495,0.709659,0.669613,0.628673,0.587177,0.545471,0.439262,0.315261,0.204999])

        #DATA

        if not self.DATASET == "mike-hires":
           raise io_mp.LikelihoodError('Error: for the time being, only the mike - hires dataset is available')
           exit()

        file_path = os.path.join(self.data_directory, self.MIKE_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           y_M_reshaped = pickle.load(pkl)
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: MIKE spectra file is missing')
           exit()
        file_path = os.path.join(self.data_directory, self.HIRES_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           y_H_reshaped = pickle.load(pkl)
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: HIRES spectra file is missing')
           exit()

        file_path = os.path.join(self.data_directory, self.MIKE_cov_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           cov_M_inverted = pickle.load(pkl)
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: MIKE covariance matrix file is missing')
           exit()

        file_path = os.path.join(self.data_directory, self.HIRES_cov_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           cov_H_inverted = pickle.load(pkl)
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: HIRES covariance matrix file is missing')
           exit()

        file_path = os.path.join(self.data_directory, self.PF_noPRACE_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           self.PF_noPRACE = pickle.load(pkl)
           pkl.close()
        else:
           raise io_mp.ConfigurationError('Error: PF_noPRACE file is missing')
           exit()

        self.cov_MH_inverted = block_diag(cov_H_inverted,cov_M_inverted)
        self.y_MH_reshaped = np.concatenate((y_H_reshaped, y_M_reshaped))

        print "Initialization of Lya likelihood done"

    #from alpha to 1./k_{1/2} in order to interpolate in a less sparse grid
    def khalf(self,alpha,beta,gamma):
        return ((((0.5)**(1/(2*gamma)) - 1)**(1/beta))/alpha)**(-1)       #1./k1half

#### functions used later on in the code, moved here to avoid nested definitions #####

    #model function #T^2=P_model/P_ref
    def T(self,k,alpha,beta,gamma):
        return (1. + (alpha*k)**(beta))**(gamma)

#    #define objective function: returns the array to be minimized # MArchi: moved fcn2min inside loglkl
#    def fcn2min(self,params, k, Tk):
#        alpha = params['alpha']
#        beta = params['beta']
#        gamma = params['gamma']
#        model = self.T(k,alpha,beta,gamma) #(1. + (alpha*kappa_interp)**(beta))**(gamma)
#        return (model - Tk)      #standard residuals

    def z_dep_func(self,parA, parS, z):  #analytical function for the redshift dependence of t0 and slope
        return parA*(( (1.+z)/(1.+self.zp) )**parS)

    def ordkrig_distance(self,p1, p2, p3, p4, p5, p6, p7, v1, v2, v3, v4, v5, v6, v7):
        return (((p1 - v1)**2 + (p2 - v2)**2 + (p3 - v3)**2 + (p4 - v4)**2 + (p5 - v5)**2 + (p6 - v6)**2 + (p7 - v7)**2)**(0.5) + self.epsilon)**self.exponent

    def ordkrig_norm(self,p1, p2, p3, p4, p5, p6, p7, v1, v2, v3, v4, v5, v6, v7):
        return np.sum(1./self.ordkrig_distance(p1, p2, p3, p4, p5, p6, p7, v1, v2, v3, v4, v5, v6, v7))

    def ordkrig_lambda(self,p1, p2, p3, p4, p5, p6, p7, v1, v2, v3, v4, v5, v6, v7):
        return (1./self.ordkrig_distance(p1, p2, p3, p4, p5, p6, p7, v1, v2, v3, v4, v5, v6, v7))/self.ordkrig_norm(p1, p2, p3, p4, p5, p6, p7, v1, v2, v3, v4, v5, v6, v7)

    def ordkrig_estimator(self,p21, z):
        pa10 = (self.z_dep_func(p21[11], p21[12], z[:])*1e4)/(self.t0_max[:]-self.t0_min[:])
        pb10 = self.z_dep_func(p21[13], p21[14], z[:])/(self.slope_max[:]-self.slope_min[:])
        p37 = np.concatenate((p21[3:11], pa10[6:], pb10[6:]))
        astrokrig_result = np.zeros((self.zeta_full_length, self.kappa_full_length ), 'float64')
        for index in range(self.num_z_XQ,len(self.redshift)):
            astrokrig_result[index,:] = np.sum(np.multiply(self.ordkrig_lambda(p37[0]/self.zreio_range, p37[1], p37[2]/self.neff_range, p37[3], p37[4+index-self.num_z_XQ]/(self.F_prior_max[index-self.num_z_overlap]-self.F_prior_min[index-self.num_z_overlap]), p37[8+index-self.num_z_XQ], p37[12+index-self.num_z_XQ], self.X[:,0], self.X[:,1], self.X[:,2], self.X[:,3], self.X[:,4+index-self.num_z_overlap], self.X[:,14+index-self.num_z_overlap], self.X[:,24+index-self.num_z_overlap]), self.input_full_matrix_interpolated_ASTRO[index,:,:]),axis=1)
        return astrokrig_result

    def ordkrig_distance_3D(self,par1, par2, par3, var1, var2, var3):
        return (((par1 - var1)**2 + (par2 - var2)**2 + (par3 - var3)**2)**(0.5) + self.epsilon)**self.exponent

    def ordkrig_norm_3D(self,par1, par2, par3):
        return np.sum(1./self.ordkrig_distance_3D(par1, par2, par3, self.X_ABG[:,0], self.X_ABG[:,1], self.X_ABG[:,2]))

    def ordkrig_lambda_3D(self,par1, par2, par3, var1, var2, var3):
        return (1./self.ordkrig_distance_3D(par1, par2, par3, var1, var2, var3))/self.ordkrig_norm_3D(par1,par2,par3)

    def ordkrig_estimator_3D(self,p21, z):
        ABG_matrix_new = np.zeros(( self.zeta_full_length, self.kappa_full_length, self.grid_size+self.num_sim_thermal), 'float64')
        NEW_ABG_matrix = np.zeros(( self.grid_size+self.num_sim_thermal, self.zeta_full_length, self.kappa_full_length), 'float64')
        full_matrix_interpolated_ABG = np.zeros(( self.zeta_full_length, self.kappa_full_length, self.grid_size+self.num_sim_thermal), 'float64')
        for i in range(self.zeta_full_length):
            for j in range(self.kappa_full_length):
                NEW_ABG_matrix[:,i,j] = self.input_full_matrix_interpolated_ABG[i,j,:]
        ABG_matrix_new = NEW_ABG_matrix + self.ordkrig_estimator(p21,z) - 1.
        ABG_matrix_new = np.clip(ABG_matrix_new, 0. , None)
        for i in range(self.zeta_full_length):
            for j in range(self.kappa_full_length):
                full_matrix_interpolated_ABG[i,j,:] = ABG_matrix_new[:,i,j]
        return np.sum(np.multiply(self.ordkrig_lambda_3D((self.khalf(p21[0],p21[1],p21[2]))/(self.a_max-self.a_min), p21[1]/(self.b_max-self.b_min), p21[2]/(self.g_max-self.g_min), self.X_ABG[:,0], self.X_ABG[:,1], self.X_ABG[:,2]), full_matrix_interpolated_ABG[:,:,:]),axis=2)

### end of block ####

    def loglkl(self, cosmo, data):

        k = np.logspace(np.log10(self.kmin), np.log10(self.kmax), num=self.k_size)

        if not os.path.exists(self.bin_file_path):
           with open(self.bin_file_path, 'w') as bin_file:
                bin_file.write('#')
                for name in data.get_mcmc_parameters(['varying']):
                    name = re.sub('[$*&]', '', name)
                    bin_file.write(' %s\t' % name)
                for name in data.get_mcmc_parameters(['derived']):
                    name = re.sub('[$*&]', '', name)
                    bin_file.write(' %s\t' % name)
                bin_file.write('\n')
                bin_file.close()

        #deal with the astro nuisance parameters
        if 'T0a' in data.mcmc_parameters:
            T0a=data.mcmc_parameters['T0a']['current']*data.mcmc_parameters['T0a']['scale']
        else:
            T0a=0.74
        if 'T0s' in data.mcmc_parameters:
            T0s=data.mcmc_parameters['T0s']['current']*data.mcmc_parameters['T0s']['scale']
        else:
            T0s=-4.38
        if 'gamma_a' in data.mcmc_parameters:
            gamma_a=data.mcmc_parameters['gamma_a']['current']*data.mcmc_parameters['gamma_a']['scale']
        else:
            gamma_a=1.45
        if 'gamma_s' in data.mcmc_parameters:
            gamma_s=data.mcmc_parameters['gamma_s']['current']*data.mcmc_parameters['gamma_s']['scale']
        else:
            gamma_s=1.93
        if 'Fz1' in data.mcmc_parameters:
            Fz1=data.mcmc_parameters['Fz1']['current']*data.mcmc_parameters['Fz1']['scale']
        else:
            Fz1=0.35
        if 'Fz2' in data.mcmc_parameters:
            Fz2=data.mcmc_parameters['Fz2']['current']*data.mcmc_parameters['Fz2']['scale']
        else:
            Fz2=0.26
        if 'Fz3' in data.mcmc_parameters:
            Fz3=data.mcmc_parameters['Fz3']['current']*data.mcmc_parameters['Fz3']['scale']
        else:
            Fz3=0.18
        if 'Fz4' in data.mcmc_parameters:
            Fz4=data.mcmc_parameters['Fz4']['current']*data.mcmc_parameters['Fz4']['scale']
        else:
            Fz4=0.07
        if 'F_UV' in data.mcmc_parameters:
            F_UV=data.mcmc_parameters['F_UV']['current']*data.mcmc_parameters['F_UV']['scale']
        else:
            F_UV=0.0

        #Temporarily store alpha beta and gamma
        if 'alpha' in self.use_nuisance:
            data.mcmc_parameters['alpha']['current']=-1.e30
        if 'beta' in self.use_nuisance:
            data.mcmc_parameters['beta']['current']=-1.e30
        if 'gamma' in self.use_nuisance:
            data.mcmc_parameters['gamma']['current']=-1.e30

        h=cosmo.h()
        Plin = np.zeros(len(k), 'float64')
        for index_k in range(len(k)):
            Plin[index_k] = cosmo.pk_lin(k[index_k]*h, 0.0)
        Plin *= h**3

        #here compute the Lya k scale
        Om=cosmo.Omega_m()
        OL=cosmo.Omega_Lambda()
        k_neff=self.k_s_over_km*100./(1.+self.z)*(((1.+self.z)**3*Om+OL)**(1./2.))
        #print 'k_neff = ',k_neff, 'h/Mpc'

        derived = cosmo.get_current_derived_parameters(data.get_mcmc_parameters(['derived']))
        for name, value in derived.iteritems():
                    data.mcmc_parameters[name]['current'] = value
        for name in derived.iterkeys():
                    data.mcmc_parameters[name]['current'] /= \
                        data.mcmc_parameters[name]['scale']

        #z_reio=cosmo.z_reio()
        z_reio=data.mcmc_parameters['z_reio']['current'] 
        if z_reio<self.zind_param_min[0]:
           z_reio = self.zind_param_min[0]
        if z_reio>self.zind_param_max[0]:
           z_reio=self.zind_param_max[0]
        sigma8=data.mcmc_parameters['sigma8']['current']
        neff=data.mcmc_parameters['neff']['current']
        #print 'z_reio = ',z_reio,'sigma8 = ',sigma8,' neff = ',neff
        #print '\n'


        #sanity check on the cosmological parameters #the one on z_reio can be removed
        if ((z_reio<self.zind_param_min[0] or z_reio>self.zind_param_max[0]) or (sigma8<self.zind_param_min[1] or sigma8>self.zind_param_max[1]) or (neff<self.zind_param_min[2] or neff>self.zind_param_max[2])):
           with open(self.bin_file_path, 'a') as bin_file:
                bin_file.write('#Error_cosmo\t')
                #for name, value in data.mcmc_parameters.iteritems():
                    #bin_file.write(' %.6e' % (value['current']*value['scale']))
                for elem in data.get_mcmc_parameters(['varying']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                for elem in data.get_mcmc_parameters(['derived']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                bin_file.write('\n')
                bin_file.close()
           sys.stderr.write('#Error_cosmo\n')
           sys.stderr.flush()
           return data.boundary_loglike

        classNeff=cosmo.Neff()

        #print '\n'
        #print 'initial model'
        #print data.cosmo_arguments
        #print '\n'

        #param_lcdm_equiv = deepcopy(data.cosmo_arguments)
        param_backup = data.cosmo_arguments.copy()

        #deal with dark radiation (ethos & Co.) according to arXiv:1412.6763
        if 'xi_idr' in data.cosmo_arguments or 'N_ur' in data.cosmo_arguments or 'N_ncdm' in data.cosmo_arguments:

            #DeltaNeff=stat_f_idr*(xi_idr**4)/7.*8./((4./11.)**(4./3.))+(N_ur+N_ncdm-3.046)
            #eta2=(1.+0.2271*(3.046+DeltaNeff))/(1.+0.2271*3.046)
            eta2=(1.+0.2271*classNeff)/(1.+0.2271*3.046)
            eta=np.sqrt(eta2)
            #print 'classNeff-3.046 = ',classNeff-3.046,' eta^2 = ',eta2
            #print '\n'

            
            if 'N_ur' in data.cosmo_arguments:
                data.cosmo_arguments['N_ur'] = 3.046 #3.0 in simulations
            if 'N_ncdm' in data.cosmo_arguments:
                del data.cosmo_arguments['N_ncdm']

            if 'omega_b' in data.cosmo_arguments:
                data.cosmo_arguments['omega_b'] *= 1./eta2
            if 'omega_cdm' in data.cosmo_arguments: #MArchi wdm might give 0*0, but this should never be the case (see eta)
                data.cosmo_arguments['omega_cdm'] *= 1./eta2
            if 'H0' in data.cosmo_arguments:
                data.cosmo_arguments['H0'] *= 1./eta
            if '100*theta_s' in data.cosmo_arguments:
                raise io_mp.ConfigurationError('Error: run with H0 instead of 100*theta_s')
                exit()

            #idm
            if 'xi_idr' in data.cosmo_arguments:
                del data.cosmo_arguments['xi_idr']
                if 'f_idm_dr' in data.cosmo_arguments:
                    if (data.cosmo_arguments['f_idm_dr'] != 1.0):
                        raise io_mp.ConfigurationError('Error: for the time being, f_idm_dr has to be set to 1.0')
                        exit()
                    del data.cosmo_arguments['f_idm_dr'] 
                #if 'a_dark' in data.cosmo_arguments:#this is not really needed, because once f_idm_dr = 0., whatever is a_dark, pba->has_idm==_FALSE_
                    #data.cosmo_arguments['a_dark'] = 0.0

            #hdm
            if 'm_ncdm' in data.cosmo_arguments and not 'omega_ncdm' in data.cosmo_arguments and not 'Omega_ncdm' in data.cosmo_arguments:
                data.cosmo_arguments['m_ncdm'] *= 1./eta2
            if 'omega_ncdm' in data.cosmo_arguments and not 'Omega_ncdm' in data.cosmo_arguments and not 'm_ncdm' in data.cosmo_arguments:
                data.cosmo_arguments['omega_ncdm'] *= 1./eta2

            #wdm
            if 'm_ncdm' in data.cosmo_arguments and ('omega_ncdm' in data.cosmo_arguments or 'Omega_ncdm' in data.cosmo_arguments):
                del data.cosmo_arguments['m_ncdm']

                del data.cosmo_arguments['k_per_decade_for_pk']
                del data.cosmo_arguments['ncdm_fluid_approximation']
                del data.cosmo_arguments['l_max_ncdm']
                del data.cosmo_arguments['Number of momentum bins']
                del data.cosmo_arguments['Maximum q']
                del data.cosmo_arguments['Quadrature strategy']

                if 'Omega_ncdm' in data.cosmo_arguments:
                    if 'Omega_cdm' in data.cosmo_arguments:
                        data.cosmo_arguments['Omega_cdm'] +=data.cosmo_arguments['Omega_ncdm']
                    if 'omega_cdm' in data.cosmo_arguments:
                        data.cosmo_arguments['omega_cdm'] +=data.cosmo_arguments['Omega_ncdm']*h*h/eta/eta #here we need the eta correction
                    del data.cosmo_arguments['Omega_ncdm']
                if 'omega_ncdm' in data.cosmo_arguments:
                    if 'Omega_cdm' in data.cosmo_arguments:
                        data.cosmo_arguments['Omega_cdm'] +=data.cosmo_arguments['omega_ncdm']/h/h #no eta correction should be applied here
                    if 'omega_cdm' in data.cosmo_arguments:
                        data.cosmo_arguments['omega_cdm'] +=data.cosmo_arguments['omega_ncdm']/eta2
                    del data.cosmo_arguments['omega_ncdm']


        cosmo.struct_cleanup()
        cosmo.empty()
        cosmo.set(data.cosmo_arguments)

        #print '\n'
        #print 'lcdm equivalent'
        #print data.cosmo_arguments
        #print '\n'

        cosmo.compute(['lensing'])

        Plin_equiv = np.zeros(len(k), 'float64')
        h = cosmo.h()
        for index_k in range(len(k)):
            Plin_equiv[index_k] = cosmo.pk_lin(k[index_k]*h, 0.0)
        Plin_equiv *= h**3

        cosmo.struct_cleanup()
        cosmo.empty()
        data.cosmo_arguments = param_backup
        cosmo.set(data.cosmo_arguments)

        #print '\n'
        #print 'back to model'
        #print data.cosmo_arguments
        #print '\n'

        cosmo.compute(['lensing'])

        Tk = np.zeros(len(k), 'float64')
        Tk = np.sqrt(Plin/Plin_equiv)

        #sanity check on the equivalent
        k_eq_der=cosmo.get_current_derived_parameters(['k_eq'])
        k_eq=k_eq_der['k_eq']/h
        #print 'k_eq',k_eq
        if any(abs(Tk[k<np.maximum(k_eq,k[0])]**2-1.0)>0.01):
        #if any(equiv_error>0.01 for x in k<k_eq):
            with open(self.bin_file_path, 'a') as bin_file:
                bin_file.write('#Error_equiv\t')
                #for name, value in data.mcmc_parameters.iteritems():
                        #bin_file.write(' %.6e' % (value['current']*value['scale']))
                for elem in data.get_mcmc_parameters(['varying']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                for elem in data.get_mcmc_parameters(['derived']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                bin_file.write('\n')
                bin_file.close()
            sys.stderr.write('#Error_equiv\n')
            sys.stderr.flush()
            return data.boundary_loglike

        spline=interpolate.splrep(k,Tk)
        der = interpolate.splev(k, spline, der=1)

        #setting k_max (i.e. cutting oscillations from the fitted region)
        for index_k in range(len(k)):
            index_k_fit_max = -1
            if Tk[index_k]<0.1 and der[index_k]>=0.: #perhaps we could find a better condition?!
               index_k_fit_max = index_k
               break
        #print k[index_k_fit_max]

        #sanity check for neff (check that the DAO do not start before k_neff)
        if k[index_k_fit_max]<k_neff:
            with open(self.bin_file_path, 'a') as bin_file:
                bin_file.write('#Error_kneff\t')
                #for name, value in data.mcmc_parameters.iteritems():
                        #bin_file.write(' %.6e' % (value['current']*value['scale']))
                for elem in data.get_mcmc_parameters(['varying']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                for elem in data.get_mcmc_parameters(['derived']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                bin_file.write('\n')
                bin_file.close()
            sys.stderr.write('#Error_kneff\n')
            sys.stderr.flush()
            return data.boundary_loglike

        #k_fit = np.zeros(len(k), 'float64')
        #Tk_fit = np.zeros(len(k), 'float64')
        k_fit = k[:index_k_fit_max]
        Tk_fit = Tk[:index_k_fit_max]

        # fitting the given linear P(k) with the {alpha,beta,gamma}-formula

        #define objective function: returns the array to be minimized
        def fcn2min(params, k, Tk):
            alpha = params['alpha']
            beta = params['beta']
            gamma = params['gamma']
            model = self.T(k,alpha,beta,gamma) #(1. + (alpha*kappa_interp)**(beta))**(gamma)
            return (model - Tk)      #standard residuals

        # create a set of Parameters
        params = Parameters()
        params.add('alpha', value=0.001, min = self.alpha_min, max = self.alpha_max) #the min and max set here are not the ones of the tables
        params.add('beta', value=2.24, min = self.beta_min, max = self.beta_max)
        params.add('gamma', value=-4.46, min= self.gamma_min, max=self.gamma_max)

        # do fit, default is with least squares method
        #t0_fit = time.clock()

        minner = Minimizer(fcn2min, params, fcn_args=(k_fit, Tk_fit))
        result = minner.minimize(method = 'leastsq')
        best_alpha = result.params['alpha'].value
        best_beta  = result.params['beta'].value
        best_gamma = result.params['gamma'].value

        #store alpha beta and gamma
        if 'alpha' in self.use_nuisance:
            data.mcmc_parameters['alpha']['current']=best_alpha/data.mcmc_parameters['alpha']['scale']
        if 'beta' in self.use_nuisance:
            data.mcmc_parameters['beta']['current']=best_beta/data.mcmc_parameters['beta']['scale']
        if 'gamma' in self.use_nuisance:
            data.mcmc_parameters['gamma']['current']=best_gamma/data.mcmc_parameters['gamma']['scale']

        #t1_fit = time.clock()

        Tk_abg=np.zeros(len(k_fit),'float64')
        Tk_abg=self.T(k_fit, best_alpha, best_beta, best_gamma)

        # write error report
        #print '\n'
        #print params
        #print '\n'
        #report_fit(result)
        #print '\n'
        #print best_alpha,best_beta,best_gamma
        #print '\n'
        #print result.chisqr, result.redchi

##        plt.xlabel('k [h/Mpc]')
##        plt.ylabel('$P_{nCDM}/P_{CDM}$')
#        plt.ylim(0.01,1.01)
#        plt.xlim(self.kmin,self.kmax)
#        plt.xscale('log')
##       plt.yscale('log')
#        plt.grid(True)
#        plt.plot(k, Tk**2, 'r')
#        plt.plot(k, (self.T(k, best_alpha, best_beta, best_gamma))**2, 'b--')
##       plt.plot(k_fit, abs(Tk_fit**2/Tk_abg**2-1.), 'k')
##       plt.show()
#        plt.savefig('grid_fit_plot.pdf')

        #sanity check on alpha beta gamma
        if ((best_alpha<self.alpha_min or best_alpha>self.alpha_max) or (best_beta<self.beta_min or best_beta>self.beta_max) or (best_gamma<self.gamma_min or best_gamma>self.gamma_max)):
           if(best_alpha<self.alpha_min or best_alpha>self.alpha_max):
               with open(self.bin_file_path, 'a') as bin_file:
                bin_file.write('#Error_a\t')
                #for name, value in data.mcmc_parameters.iteritems():
                    #bin_file.write(' %.6e' % (value['current']*value['scale']))
                for elem in data.get_mcmc_parameters(['varying']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                for elem in data.get_mcmc_parameters(['derived']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                bin_file.write('\n')
                bin_file.close()
               sys.stderr.write('#Error_a\n')
               sys.stderr.flush()
           else:
               with open(self.bin_file_path, 'a') as bin_file:
                bin_file.write('#Error_bg\t')
                #for name, value in data.mcmc_parameters.iteritems():
                    #bin_file.write(' %.6e' % (value['current']*value['scale']))
                for elem in data.get_mcmc_parameters(['varying']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                for elem in data.get_mcmc_parameters(['derived']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                bin_file.write('\n')
                bin_file.close()
               sys.stderr.write('#Error_bg\n')
               sys.stderr.flush()

           return data.boundary_loglike

        #sanity check on the alpha beta gamma fit wrt the model
        #fit_error=abs(Tk_fit/Tk_abg-1.)
        #if any(x>0.1 for x in fit_error):

        for ik in range(len(k_fit)):
            index_k_check_max = -1
            if Tk_fit[ik]<0.2:
               index_k_check_max = ik
               break

        if any(abs(Tk_fit[:index_k_check_max]/Tk_abg[:index_k_check_max]-1.)>0.1):
            with open(self.bin_file_path, 'a') as bin_file:
                bin_file.write('#Error_fit\t')
                #for name, value in data.mcmc_parameters.iteritems():
                     #bin_file.write(' %.6e' % (value['current']*value['scale']))
                for elem in data.get_mcmc_parameters(['varying']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                for elem in data.get_mcmc_parameters(['derived']):
                    bin_file.write(' %.6e\t' % data.mcmc_parameters[elem]['current'])
                bin_file.write('\n')
                bin_file.close()
            sys.stderr.write('#Error_fit\n')
            sys.stderr.flush()
            return data.boundary_loglike



        chi2=0.
        #theta=np.zeros(len(self.use_nuisance)+self.params_numbers+len(self.zind_param_size)-1, 'float64')
        model_H = np.zeros (( len(self.zeta_range_mh), len(self.k_mh) ), 'float64')
        #y_H = np.zeros (( len(self.zeta_range_mh), len(self.k_mh) ), 'float64')
        model_M = np.zeros (( len(self.zeta_range_mh)-1, len(self.k_mh) ), 'float64')
        #y_M = np.zeros (( len(self.zeta_range_mh)-1, len(self.k_mh) ), 'float64')
        theta=np.array([best_alpha,best_beta,best_gamma,z_reio,sigma8,neff,F_UV,Fz1,Fz2,Fz3,Fz4,T0a,T0s,gamma_a,gamma_s])
        model = self.PF_noPRACE*self.ordkrig_estimator_3D(theta, self.redshift_list)
        upper_block = np.vsplit(model, [7,11])[0]
        lower_block = np.vsplit(model, [7,11])[1]

        model_H[:,:] = lower_block[:,19:]
        model_H_reshaped = np.reshape(model_H, -1, order='C')
        model_M[:,:] = lower_block[:3,19:]
        model_M_reshaped = np.reshape(model_M, -1, order='C')
        model_MH_reshaped = np.concatenate((model_H_reshaped,model_M_reshaped))
        chi2 = np.dot((self.y_MH_reshaped - model_MH_reshaped),np.dot(self.cov_MH_inverted,(self.y_MH_reshaped - model_MH_reshaped)))

        loglkl = - 0.5 * chi2

        #print 'Lya -chi^2/2 = ',loglkl

        return loglkl
