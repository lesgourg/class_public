from genericpath import exists
import multiprocessing
import random
import copy
import numpy as np
from scipy.interpolate import CubicSpline
import os
import h5py as h5

from tqdm import tqdm
import classy
from classy import Class

from classynet.generate.generate_cosmological_parameters import sample_cosmological_parameters
from classynet.plotting.plot_tester_cls import *
import classynet.utils as utils


'''
The Tester class allows to obtain and plot all kinds of data and information from the ClassNet
'''
class Tester:
    def __init__(self, workspace):
        self.workspace = workspace

        self.k = self.workspace.loader().k()
        manifest = workspace.loader().manifest()
        self.fixed = manifest["fixed"]

        self.WDIR = workspace.path


    
    def create_cls(self,pnames,N=1000):
        """
        This function generates spectra with and without the neural network of the test set. 
        It will be used to qualify the differences between ClassFULL and ClassNet.
        The networks have to be trained for this purpose
        The test/parameters.h5 hasn't been sampled in advance it is going to be here.
        """
        if not os.path.isfile(self.workspace.test_data / 'cl_test_parameters.h5'):
            domain = self.workspace.domain_from_path(pnames)

            #Sample parameters according to the domain and save them in the workspace as "samples.h5"
            domain.sample_save(training_count=0, validation_count=0, test_count=N, file_name='cl_test_parameters')

        # Load parameters
        training, validation, test = self.workspace.loader().cosmological_parameters(file_name='cl_test_parameters')
        #Create the spectra data for class_full
        if not os.path.isfile(self.workspace.test_data / 'FULL_cls.h5'):
            self.workspace.generator().generate_spectra_data(
                #training=training,
                #validation=validation,
                test=test,
                processes=4,
                use_nn=False,
            )

        #Create the spectra data for class_full
        if not os.path.isfile(self.workspace.test_data / 'NN_cls.h5'):
            self.workspace.generator().generate_spectra_data(
                #training=training,
                #validation=validation,
                test=test,
                processes=4,
                use_nn=True,
            )

    def load_files(self,path,my_set):
        with h5.File(path, "r") as f:
            return {key: list(f[my_set].get(key)) for key in f[my_set].keys()}

    def compare_cls_full_net(self,pnames=None,mode='cl_lens',N_lines = 100, N_contain = 1000):
        """
        This function loades the cls stored in '/test/FULL_cls.h5' and '/test/NN_cls.h5' and compares the differences and make plots.
        """

        # load cls and parameters
        FULL_cls = self.load_files(self.workspace.test_data / 'FULL_cls.h5',mode)
        NN_cls = self.load_files(self.workspace.test_data / 'NN_cls.h5',mode)
        NN_para = self.load_files(self.workspace.test_data / 'NN_cls.h5','parameter')
        FULL_para = self.load_files(self.workspace.test_data / 'FULL_cls.h5','parameter')

        # get the chisquares of the individual data points
        if pnames != None:
            domain = self.workspace.domain_from_path(pnames)
        else:
            domain = self.workspace.domain()

        FULL_para_list_of_dirs = utils.transpose_dict_of_lists(FULL_para)

        chisq_list = []
        for i,dir in enumerate(FULL_para_list_of_dirs):
            _, chisq = domain.contains(dir)
            chisq_list.append(chisq)

        chisq_list = np.array(chisq_list)

        # spectra which should be plotted
        spectra = ['tt','ee','bb','te','tp','pp','pk']

        # Start with line plots which show ClassNN and CLassFULL in 1 plot
        os.makedirs(self.workspace.plots / 'cl_tests' / 'line_plots', exist_ok=True)
        os.makedirs(self.workspace.plots / 'cl_tests' / 'variance_plots', exist_ok=True)
        os.makedirs(self.workspace.plots / 'cl_tests' / 'contain_plots', exist_ok=True)

        my_plotter = CL_plotter()
        for spectrum in spectra:
            print('plotting {}'.format(spectrum))
            save_dir = self.workspace.plots / 'cl_tests' / 'line_plots' / (spectrum + '.pdf')
            my_plotter.line_plots(FULL_cls,NN_cls,spectrum,save_dir,N=N_lines)
            save_dir = self.workspace.plots / 'cl_tests' / 'line_plots' / (spectrum + '_dif_abs.pdf')
            my_plotter.line_plots_difference(FULL_cls,NN_cls,spectrum,save_dir,measure="abs",N=N_lines)
            save_dir = self.workspace.plots / 'cl_tests' / 'line_plots' / (spectrum + '_dif_abs_chi.pdf')
            my_plotter.line_plots_difference(FULL_cls,NN_cls,spectrum,save_dir,measure="abs",chisq = chisq_list,N=N_lines)
            save_dir = self.workspace.plots / 'cl_tests' / 'line_plots' / (spectrum + '_dif_rel.pdf')
            my_plotter.line_plots_difference(FULL_cls,NN_cls,spectrum,save_dir,measure="rel",N=N_lines)


            save_dir = self.workspace.plots / 'cl_tests' / 'variance_plots' / (spectrum + '_cosmic_variance.pdf')
            my_plotter.line_plots(FULL_cls,NN_cls,spectrum,save_dir,cosmic_variance=True,N=N_lines)

            sigmas=[0.99,0.95,0.68]
            save_dir = self.workspace.plots / 'cl_tests' / 'contain_plots' / (spectrum + '_dif_cont_error.pdf')
            my_plotter.line_plots_contain(FULL_cls,NN_cls,spectrum,save_dir,noise=True,sigmas = sigmas, N=N_contain)
            save_dir = self.workspace.plots / 'cl_tests' / 'contain_plots' / (spectrum + '_dif_cont.pdf')
            my_plotter.line_plots_contain(FULL_cls,NN_cls,spectrum,save_dir,noise=False,sigmas = sigmas, N=N_contain)



























    def test(self, count=None, cheat=None, seed=None,nonlinear=None):
        """
        test the performance of the trained networks by computing the observables
        (Cl's and P(k)) for `count` cosmological scenarios.
        For each observable a plot will be created that displays the deviation
        of the network prediction from the true value of the observable.

        Optionally, `cheat` can be provided as a list of names of source functions
        for which to ignore the network predictions and instead take the true
        value from CLASS.
        This aids in isolating networks which may not perform well.
        """

        _, params = self.workspace.loader().cosmological_parameters()
        validation_size = len(params[next(iter(params))])
        if count is None:
            count = validation_size

        if count > validation_size:
            print("WARNING: You requested a validation same count of {0}, \
                  but there are only {1}! Proceeding with {1}".format(count, validation_size))
            count = validation_size

        if seed is not None:
            random.seed(seed)
        selection = random.sample(range(validation_size), count)
        print("Evaluating network predictions for {} cosmologies.".format(count))

        cosmo_params = [{k: v[i] for k, v in params.items()} for i in selection]

        # TODO TODO TODO REMOVE THIS LINE!!!!!
        # cosmo_params = sorted(cosmo_params, key=lambda item: np.abs(item["Omega_k"] - (-0.011))) # removed this line, but unsure if more is required to compensate
        # TODO TODO TODO REMOVE THIS LINE!!!!!

        # params = sample_cosmological_parameters(self.domain, count)

        # class_params = list(map(self.get_class_parameters, cosmo_params))
        # class_params_nn = list(map(self.get_class_parameters_nn, cosmo_params))
        class_params    = [self.get_class_parameters(p,nonlinear=nonlinear) for p in cosmo_params]
        class_params_nn = [self.get_class_parameters_nn(p, cheat=cheat,nonlinear=nonlinear) for p in cosmo_params]
        print("class_params",class_params)
        print("class_params_nn",class_params_nn)
        # for recording cl/pk errors as functions of cosmological parameters
        stats = []

        class_pair_iter = map(self.run_class_for, zip(class_params, class_params_nn))
        counter = 0
        for cosmo_param_dict, (exc, pair) in zip(cosmo_params, class_pair_iter):
            print("Omega_k", cosmo_param_dict["Omega_k"])
            if exc:
                print("An exception occured:")
                print(str(exc))
                print("continuing...")
                continue
            (cl_true, k_pk_true, pk_true, k, dm_true, cl_lens_true), (cl_nn, k_pk_nn, pk_nn, k_nn, dm_nn, cl_lens_nn) = pair
            print("cltrue-clnn",cl_true["tt"]-cl_nn["tt"])
            counter += 1
            print("PROGRESS: {}".format(counter))

            self.update_stats(stats, cosmo_param_dict, cl_true, k, dm_true, cl_nn, k_pk_true, pk_true, k_pk_nn, pk_nn, k_nn, dm_nn, cl_lens_true, cl_lens_nn)

        self.save_stats(stats)

    def update_stats(self, stats, cosmo_params, cl_true, k, dm_true, cl_nn, k_pk_true, pk_true, k_pk_nn, pk_nn, k_nn, dm_nn,cl_lens_true,cl_lens_nn):
        cl_err = {key: cl_nn[key] - cl_true[key] for key in cl_true}

        # since P_NN(k) and P_true(k) _may_ be sampled on different k grids, we
        # need to interpolate (in this case, onto the k_pk_true)
        pk_spline = CubicSpline(k_pk_nn, pk_nn)
        pk_nn_resampled = pk_spline(k_pk_true)
        pk_err = pk_nn_resampled - pk_true

        cl_err_relative = {key: (cl_nn[key] - cl_true[key]) / cl_true[key] for key in cl_true}
        pk_err_relative = pk_err / pk_true

        stat_dict = {
            "parameters": cosmo_params,
            "cl_true": cl_true,
            "cl_nn": cl_nn,
            "cl_error": cl_err,
            "k_pk": k_pk_true,
            "pk_true": pk_true,
            "k_pk_nn": k_pk_nn,
            "pk_nn": pk_nn,
            "pk_error": pk_err,
            "cl_error_relative": cl_err_relative,
            "pk_error_relative": pk_err_relative,
            "k": k,
            "k_nn": k_nn,
            "delta_m": dm_true,
            "delta_m_nn": dm_nn,
            "cl_lens_true":cl_lens_true,
            "cl_lens_nn":cl_lens_nn,
        }

        # stat_dict = {k: v if not isinstance(v, np.ndarray) else list(v) for k, v in stat_dict.items()}
        stats.append(stat_dict)

    def get_class_parameters(self, cosmo_params,nonlinear=None):
        """ add fixed parameters (no nn) """
        params = {}
        params.update(cosmo_params)
        params.update(self.fixed)
        if nonlinear!= None and "non linear" in params:
            del params["non linear"]
        if nonlinear=="halofit":
            params.update({"non_linear":"halofit"})
        elif nonlinear=="linear" and "non_linear" in params:
            del params["non_linear"]
        return params

    def get_class_parameters_nn(self, cosmo_params, cheat=None,nonlinear=None):
        """ add fixed parameters (with nn) """
        params = self.get_class_parameters(cosmo_params)
        params["neural network path"] = str(self.workspace.path.resolve())
        params["nn_verbose"] = 1
        if nonlinear!= None and "non linear" in params:
            del params["non linear"]
        if nonlinear=="halofit":
            params.update({"non_linear":"halofit"})
        elif nonlinear=="linear" and "non_linear" in params:
            del params["non_linear"]
        # if cheating is enabled, we need to notify classy of the quantities
        if cheat:
            params["nn_cheat"] = cheat
        return params

    def parameter_dicts(self, cosmo_params, cheat=None):
        """
        iterator over tuples of `(params_without_nn, params_with_nn)` for
        given domain.
        """
        count = len(cosmo_params)
        for i in range(count):
            params = {}
            params.update(self.fixed)
            params.update({k: v[i] for k, v in cosmo_params.items()})

            params_nn = params.copy()
            params_nn["neural network path"] = str(self.workspace.path.resolve())
            # if cheating is enabled, we need to notify classy of the quantities
            if cheat:
                params_nn["nn_cheat"] = cheat

            yield params, params_nn

    def run_class_for(self, param_tuple):
        """
        run class twice: once without nn, and once with.
        """
        params, params_nn = param_tuple
        try:
            ret = self.run_class(params), self.run_class(params_nn)
            exc = None
        except (classy.CosmoSevereError, classy.CosmoComputationError) as e:
            ret = None
            exc = e

        return exc, ret

    def k_pk(self, cosmo):
        k_min = cosmo.k_min()
        k = self.k[self.k > k_min]
        # TODO don't hardcode limit?
        k = np.insert(k, 0, k_min)
        k = k[k < 100]
        return k

    def run_class(self, params):
        """
        Run CLASS for the given `params` and return the Cls and mPk.
        """
        import time
        cosmo = Class()
        cosmo.set(params)
        print(cosmo.pars)
        start = time.perf_counter()
        report = {}
        cosmo.compute(performance_report=report)
        end = time.perf_counter()
        elapsed = end - start
        print("running class took {}s".format(elapsed))
        cls = truncate(cosmo.raw_cl())
        print("l raw:\n",cls["ell"])
        cls_lens = truncate(cosmo.lensed_cl())
        print("l lensed:\n",cls_lens["ell"])
        k_pk = self.k_pk(cosmo)
        # TODO remove assertion in release
        assert np.all(k_pk >= cosmo.k_min())
        pk = np.vectorize(cosmo.pk)(k_pk, 0.0)
        
        sources, k, tau = cosmo.get_sources()

        # also take delta_m @ today
        k = k.copy()
        delta_m = sources["delta_m"][:, -1].copy()
        ret = copy.deepcopy(cls), copy.deepcopy(k_pk), copy.deepcopy(pk), copy.deepcopy(k), copy.deepcopy(delta_m), copy.deepcopy(cls_lens)

        cosmo.struct_cleanup()

        return ret

    def save_stats(self, stats):
        import pickle
        path = self.workspace.stats_file
        print("Writing Cl/Pk error stats to '{}'".format(path))
        with open(path, "wb") as out:
            pickle.dump(stats, out)









def truncate(cls):
    """
    Removes the first two multipoles of a Cl dict.
    """
    return {k: v[2:] for k, v in cls.items()}



