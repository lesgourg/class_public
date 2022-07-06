"""
.. module:: classy
    :synopsis: Python wrapper around CLASS
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>

This module defines a class called Class. It is used with Monte Python to
extract cosmological parameters.

# JL 14.06.2017: TODO: check whether we should free somewhere the allocated fc.filename and titles, data (4 times)

"""
import time
import os
from math import exp,log
import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
from cython.parallel import prange
cimport cython
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from cython.view cimport array as cvarray

import classynet.workspace
import classynet.predictors
import classynet.models

import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()
def string_to_char(d):
    if (sys.version_info >= (3,0) and isinstance(d,str)) or isinstance(d,unicode):
        d = d.encode('utf8')
    return d

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i



# Import the .pxd containing definitions
from cclassy cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

__version__ = _VERSION_.decode("utf-8")

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass

class Timer:
    """
    Simple help for performance measurements.
    """
    def __init__(self):
        self._start = {}
        self._end = {}
        self._times = {}

    def start(self, name):
        if name in self._start:
            print("WARNING: Overwriting measurement {}".format(name))
        self._start[name] = time.perf_counter()

    def end(self, name):
        if name not in self._start:
            raise ValueError(
               "Measurement '{}' has not started; cannot end!".format(name)
               )
        if name in self._end:
            print("WARNING: Overwriting measurement {}".format(name))
        self._end[name] = time.perf_counter()
        self._times[name] = self._end[name] - self._start[name]

    @property
    def times(self):
        return self._times

    def __getitem__(self, name):
        return self._times[name]

    def __setitem__(self, name, value):
        self._times[name] = value

cdef class Class:
    """
    Class wrapping, creates the glue between C and python

    The actual Class wrapping, the only class we will call from MontePython
    (indeed the only one we will import, with the command:
    from classy import Class

    """
    # List of used structures, defined in the header file. They have to be
    # "cdefined", because they correspond to C structures
    cdef precision pr
    cdef background ba
    cdef thermodynamics th
    cdef perturbations pt
    cdef primordial pm
    cdef fourier fo
    cdef transfer tr
    cdef harmonic hr
    cdef output op
    cdef lensing le
    cdef distortions sd
    cdef file_content fc
    #cdef np.ndarray[np.double_t, ndim=2, mode="c"] NN_prediction_numpy #SG cannot defined here

    cpdef int computed # Flag to see if classy has already computed with the given pars
    cpdef int allocated # Flag to see if classy structs are allocated already
    cpdef object _pars # Dictionary of the parameters
    cpdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

    cpdef bint using_NN
    cpdef bint pred_init
    cpdef dict NN_generations
    cpdef object predictor
    cpdef list NN_source_index_list
    cpdef list NN_pnames


    cdef double [:, :] NN_prediction # The NN output which will be handed over to CLASS
    cpdef int NN_source_size

    # Defining two new properties to recover, respectively, the parameters used
    # or the age (set after computation). Follow this syntax if you want to
    # access other quantities. Alternatively, you can also define a method, and
    # call it (see _T_cmb method, at the very bottom).
    property pars:
        def __get__(self):
            return self._pars
    property state:
        def __get__(self):
            return True
    property Omega_nu:
        def __get__(self):
            return self.ba.Omega0_ncdm_tot
    property nonlinear_method:
        def __get__(self):
            return self.fo.method

    def set_default(self):
        _pars = {
            "output":"tCl mPk",}
        self.set(**_pars)

    def __cinit__(self, default=False):
        cpdef char* dumc
        self.allocated = False

        # NN stuff
        self.using_NN = False
        self.NN_generations = {}
        self.pred_init = False
        self.NN_source_index_list = []
        self.NN_pnames = []

        self.computed = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()

    def __dealloc__(self):
        if self.allocated:
          self.struct_cleanup()
        self.empty()
        # Reset all the fc to zero if its not already done
        if self.fc.size !=0:
            self.fc.size=0
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
            free(self.fc.filename)

    # Set up the dictionary
    def set(self,*pars,**kars):
        oldpars = self._pars.copy()
        if len(pars)==1:
            self._pars.update(dict(pars[0]))
        elif len(pars)!=0:
            raise CosmoSevereError("bad call")
        self._pars.update(kars)
        if viewdictitems(self._pars) <= viewdictitems(oldpars):
            return # Don't change the computed states, if the new dict was already contained in the previous dict
        self.computed=False
        return True

    def use_nn(self):
        """
        Utility methods that returns whether neural networks are enabled
        by checking whether 'workspace_path' is in the input parameters.
        Also checks the 'nn_verbose' parameter to determine how much information 
        to print about the usage of neural networks.
        """

        if "nn_verbose" in self.pars:
            nn_verbose = self.pars["nn_verbose"]
        else:
            nn_verbose = 0
        if not "use_nn" in self.pars:
            return False
        else:
            if not self.pars["use_nn"].lower()=='yes':
                return False
        
        if not "workspace_path" in self._pars:
                raise ValueError("'use_nn' is yes but no workspace_path is not provided!")
        elif "workspace_path" in self.pars:
            self.using_NN = self.can_use_nn()
            if nn_verbose>1:
                if not self.using_NN:
                    print("##################################")
                    print("#   NOT USING NEURAL NETWORKS!   #")
                    print("##################################")
                    return False
                else:
                    print("##################################")
                    print("#    USING NEURAL NETWORKS!      #")
                    print("##################################")
                    return True
            elif nn_verbose==1:
                if not self.using_NN:
                    print("USING NEURAL NETWORKS : False!")
                    return False
                else:
                    print("USING NEURAL NETWORKS : True!")
                    return True
            elif nn_verbose==0:
                if not self.using_NN:
                    return False
                else:
                    return True
            else:
                if not self.using_NN:
                    print("##################################")
                    print("#   NOT USING NEURAL NETWORKS!   #")
                    print("##################################")
                    return False
                else:
                    print("##################################")
                    print("#    USING NEURAL NETWORKS!      #")
                    print("##################################")
                    return True


    # This function checks for the cosmological parameters whether classnet can be used
    def can_use_nn(self):
        """ may only be called if neural networks are enabled """        
        workspace = self.nn_workspace()
        domain = workspace.loader().domain_descriptor()

        # load the relevant cosmological parameter and check whether they lay indide the domain
        nn_cosmo_pars = self.nn_cosmological_parameters()
        domain_use_nn, self.pt.network_deltachisquared = domain.contains(nn_cosmo_pars)

        # if 'nn_force' was set to 'yes' we will overwrite the recommendation and force classnet to use the neural networks
        if 'nn_force' in self._pars:
            if self._pars['nn_force'].lower() == 'yes':
                domain_use_nn=True

        if not domain_use_nn:
            if "nn_verbose" in self.pars:
                if self.pars["nn_verbose"]>1:
                    print("neural network domain of validity does not contain requested parameters due to the delta chisquare of "+ str(self.pt.network_deltachisquared))
            return False

        if self.ba.N_ncdm!=1:
            raise CosmoSevereError("You use classnet, but did not set 'N_ncdm':1")
        if self.ba.deg_ncdm[0]!=3:
            raise CosmoSevereError("You use classnet, but did not set 'deg_ncdm':3")
        if self.ba.Omega0_lambda!=0:
            raise CosmoSevereError("You use classnet, but did not set 'Omega_Lambda':0")
        if self.th.compute_damping_scale!=1:
            raise CosmoSevereError("You use classnet, but did not set 'compute damping scale':'yes'")
        pk_max = nn_cosmo_pars.get("P_k_max_1/Mpc")
        if pk_max is not None and pk_max > 100.0:
            raise CosmoSevereError("neural networks only applicable with 'P_k_max_1/Mpc' <= 100.0")

        return True

    def empty(self):
        self._pars = {}
        self.computed = False

    # this function collects the relevant cosmological parameters for the classnet evaluation
    def nn_cosmological_parameters(self):
        output = {}
        output['omega_b'] = self.ba.Omega0_b * self.ba.h**2
        output['omega_cdm'] = self.ba.Omega0_cdm * self.ba.h**2
        output['h'] = self.ba.h
        output['w0_fld'] = self.ba.w0_fld
        output['wa_fld'] = self.ba.wa_fld
        output['omega_ncdm'] = self.ba.Omega0_ncdm_tot * self.ba.h**2
        output['Omega_k'] = self.ba.Omega0_k
        output['tau_reio'] = self.th.tau_reio
        output['N_ur'] = ( self.ba.Neff - 3.044 ) + 0.00641
        return(output)

    # Create an equivalent of the parameter file. Non specified values will be
    # taken at their default (in Class)
    def _fillparfile(self):
        cdef char* dumc

        if self.fc.size!=0:
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
        self.fc.size = len(self._pars)
        self.fc.name = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.name!=NULL)

        self.fc.value = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.value!=NULL)

        self.fc.read = <short*> malloc(sizeof(short)*len(self._pars))
        assert(self.fc.read!=NULL)

        # fill parameter file
        i = 0
        for kk in self._pars:

            dumcp = kk.encode()
            dumc = dumcp
            sprintf(self.fc.name[i],"%s",dumc)
            dumcp = str(self._pars[kk]).strip().encode()
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = _FALSE_
            i+=1

    # Called at the end of a run, to free memory 
    def struct_cleanup(self):
        if(self.allocated != True):
          return
        if self.using_NN:
            # allocate the source functions which were overwritten/freed by the NN. They are to be freed by pert_free()
            for indices in self.NN_source_index_list:
                self.pt.sources[indices[0]][indices[1]] = <double *> malloc(sizeof(double))
            self.NN_source_index_list = []
            self.using_NN = False
        if "distortions" in self.ncp:
            distortions_free(&self.sd)
        if "lensing" in self.ncp:
            lensing_free(&self.le)
        if "harmonic" in self.ncp:
            harmonic_free(&self.hr)
        if "transfer" in self.ncp:
            transfer_free(&self.tr)
        if "fourier" in self.ncp:
            fourier_free(&self.fo)
        if "primordial" in self.ncp:
            primordial_free(&self.pm)
        if "perturb" in self.ncp:
            perturbations_free(&self.pt)
        if "thermodynamics" in self.ncp:
            thermodynamics_free(&self.th)
        if "background" in self.ncp:
            background_free(&self.ba)        
        self.allocated = False
        self.computed = False
        
    def _check_task_dependency(self, level):
        """
        Fill the level list with all the needed modules

        .. warning::

            the ordering of modules is obviously dependent on CLASS module order
            in the main.c file. This has to be updated in case of a change to
            this file.

        Parameters
        ----------

        level : list
            list of strings, containing initially only the last module required.
            For instance, to recover all the modules, the input should be
            ['lensing']

        """
        if "distortions" in level:
            if "lensing" not in level:
                level.append("lensing")
        if "lensing" in level:
            if "harmonic" not in level:
                level.append("harmonic")
        if "harmonic" in level:
            if "transfer" not in level:
                level.append("transfer")
        if "transfer" in level:
            if "fourier" not in level:
                level.append("fourier")
        if "fourier" in level:
            if "primordial" not in level:
                level.append("primordial")
        if "primordial" in level:
            if "perturb" not in level:
                level.append("perturb")
        if "perturb" in level:
            if "thermodynamics" not in level:
                level.append("thermodynamics")
        if "thermodynamics" in level:
            if "background" not in level:
                level.append("background")
        if len(level)!=0 :
            if "input" not in level:
                level.append("input")
        return level

    def _pars_check(self, key, value, contains=False, add=""):
        val = ""
        if key in self._pars:
            val = self._pars[key]
            if contains:
                if value in val:
                    return True
            else:
                if value==val:
                    return True
        if add:
            sep = " "
            if isinstance(add,str):
                sep = add

            if contains and val:
                    self.set({key:val+sep+value})
            else:
                self.set({key:value})
            return True
        return False

    cdef void overwrite_source_function(self, int index_md, int index_ic,
            int index_type,
            int k_NN_size, int tau_size, 
            double[:] S):
        """
        This utility function overwrites a single source function specified by `index_type`
        with the given source function `S`.
        Used by NN code (see below in "perturb" section of `compute()`).
        """
        cdef int tp_size = self.pt.tp_size[index_md]
        cdef int index_tau
        cdef int index_k
        ## -> Not reuqired if perform_NN_skip is true :
        # Required again because all source functions have been allocated earlier
        free(self.pt.sources[index_md][index_ic * tp_size + index_type])
        
        # Hand over the address of the numpy data array
        self.pt.sources[index_md][index_ic * tp_size + index_type] = &S[0]

        # Add the overwritten indices to list such that it can be reallocated when they are cleaned up by pert_free()
        self.NN_source_index_list.extend([[index_md, index_ic * tp_size + index_type]])

    def compute(self, level=["distortions"], performance_report=None, post_perturb_callback=None):
        """
        compute(level=["distortions"])

        Main function, execute all the _init methods for all desired modules.
        This is called in MontePython, and this ensures that the Class instance
        of this class contains all the relevant quantities. Then, one can deduce
        Pk, Cl, etc...

        Parameters
        ----------
        level : list
                list of the last module desired. The internal function
                _check_task_dependency will then add to this list all the
                necessary modules to compute in order to initialize this last
                one. The default last module is "lensing".

        .. warning::

            level default value should be left as an array (it was creating
            problem when casting as a set later on, in _check_task_dependency)

        """
        cdef ErrorMsg errmsg
        cdef int i

        timer = Timer()

        timer.start("compute")

        # Append to the list level all the modules necessary to compute.
        level = self._check_task_dependency(level)

        # Check if this function ran before (self.computed should be true), and
        # if no other modules were requested, i.e. if self.ncp contains (or is
        # equivalent to) level. If it is the case, simply stop the execution of
        # the function.
        if self.computed and self.ncp.issuperset(level):
            return

        # Check if already allocated to prevent memory leaks
        if self.allocated:
            self.struct_cleanup()

        # Otherwise, proceed with the normal computation.
        self.computed = False

        # Equivalent of writing a parameter file
        self._fillparfile()

        # self.ncp will contain the list of computed modules (under the form of
        # a set, instead of a python list)
        self.ncp=set()
        # Up until the empty set, all modules are allocated
        # (And then we successively keep track of the ones we allocate additionally)
        self.allocated = True

        # --------------------------------------------------------------------
        # Check the presence for all CLASS modules in the list 'level'. If a
        # module is found in level, executure its "_init" method.
        # --------------------------------------------------------------------
        # The input module should raise a CosmoSevereError, because
        # non-understood parameters asked to the wrapper is a problematic
        # situation.
        if "input" in level:
            timer.start("input")
            if input_read_from_file(&self.fc, &self.pr, &self.ba, &self.th,
                                    &self.pt, &self.tr, &self.pm, &self.hr,
                                    &self.fo, &self.le, &self.sd, &self.op, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.ncp.add("input")
            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            # GS: added this because neural network arguments are not relevant
            # to the C code.
            problematic_exceptions = set(["workspace_path", "nn_force", "nn_verbose","use_nn","Net_phi_plus_psi","Net_ST0_ISW","Net_ST0_Reco","Net_ST0_Reio","Net_ST1","Net_ST2_Reco","Net_ST2_Reio"])
            for i in range(self.fc.size):
                if self.fc.read[i] == _FALSE_:
                    name = self.fc.name[i].decode()
                    # GS: if parameter is an exception, do not raise problem flag
                    if name in problematic_exceptions:
                        continue
                    problem_flag = True
                    problematic_parameters.append(name)
            if problem_flag:
                raise CosmoSevereError(
                    "Class did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))
            timer.end("input")

        # The following list of computation is straightforward. If the "_init"
        # methods fail, call `struct_cleanup` and raise a CosmoComputationError
        # with the error message from the faulty module of CLASS.
        if "background" in level:
            timer.start("background")
            if background_init(&(self.pr), &(self.ba)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.ba.error_message)
            self.ncp.add("background")
            timer.end("background")

        if "thermodynamics" in level:
            timer.start("thermodynamics")
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.th.error_message)
            self.ncp.add("thermodynamics")
            timer.end("thermodynamics")

        # define objects for NN
        cdef:
            int i_index_type, index_k, index_tau, i_k
            int k_NN_size, tau_NN_size
            double [:] k_CLASS
            double [:] tau_CLASS
            # double * c_tau_NN
            # int c_tau_NN_size
            # double [:] numpy_k_CLASS
            # double [:] numpy_tau_CLASS
            int index_md
            int k_size
            int tau_size
            int index_ic
            int index_type
            # int tp_size = 0
            int tp_size
            int tot_num_of_sources = 11
            int index_tp_x
            # double [:,:] NN_interpolated
            # double [:] NN_interpolated
            # TODO remove some of the unused ones here
            #double [:, :] NN_prediction
            double tau1=0.0
            #np.ndarray[np.double_t, ndim=2, mode="c"] NN_prediction_numpy

        if "perturb" in level:

            timer.start("perturb")
            timer.start("perturb_init")


            # Allocate memory for ALL source functions (since transfer.c iterates over them)
            self.using_NN = self.use_nn()   #0.0007 sec

            if self.using_NN:
                if "nn_verbose" in self.pars:
                    if self.pars["nn_verbose"]>2:
                        print("Using neural networks; skipping regular perturbation module.")
                self.pt.perform_NN_skip = _TRUE_

            start = time.time()
            if perturbations_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            self.ncp.add("perturb")

            timer.end("perturb_init")

            # flag for using NN
            if self.using_NN:
                timer.start("neural network complete")

                timer.start("update predictor")
                # Check whether the predictor has been already build, otherwise build it now. If it already has been build, the cosmo parameter have to be updated!
                if self.pred_init==False:
                    self.predictor = classynet.predictors.build_predictor(self)
                    self.predictor.update_predictor(self)
                    self.pred_init = True
                else:
                    self.predictor.update_predictor(self)
                timer.end("update predictor")   # 1e-5 sec

                timer.start("neural network initialization")

                # Obtain prediction parameters
                index_md = self.pt.index_md_scalars;
                k_size = self.pt.k_size[index_md];
                tau_size = self.pt.tau_size;
                index_ic = self.pt.index_ic_ad;

                tp_size = self.pt.tp_size[index_md];

                k_NN = self.predictor.get_k()
                k_NN_size = len(k_NN)

                #load requested tau values
                tau_CLASS = np.zeros((tau_size))
                for index_tau in range(tau_size):
                    tau_CLASS[index_tau] = self.pt.tau_sampling[index_tau]

                NN_requested_index_types = []
                
                # add all sources that class calculates to the list for predicting
                source_names = []
                if self.pt.has_source_t:
                    NN_requested_index_types.extend([
                        self.pt.index_tp_t0,
                        self.pt.index_tp_t1,
                        self.pt.index_tp_t2,
                        ])
                    source_names.extend(["t0", "t1", "t2"])

                if self.pt.has_source_phi_plus_psi:
                    NN_requested_index_types.append(self.pt.index_tp_phi_plus_psi)
                    source_names.append("phi_plus_psi")

                if self.pt.has_source_delta_m:
                    NN_requested_index_types.append(self.pt.index_tp_delta_m)
                    source_names.append("delta_m")

                if self.pt.has_source_delta_cb:
                    NN_requested_index_types.append(self.pt.index_tp_delta_cb)
                    source_names.append("delta_cb")

                if self.pt.has_source_p:
                    NN_requested_index_types.append(self.pt.index_tp_p)
                    source_names.append("t2_p")
                    
                if self.nn_debug_enabled():
                    NN_requested_index_types.extend([
                        self.pt.index_tp_t0_reco_no_isw,
                        self.pt.index_tp_t0_reio_no_isw,
                        self.pt.index_tp_t0_isw,
                        self.pt.index_tp_t2_reco,
                        self.pt.index_tp_t2_reio,
                        ])
                    source_names.extend([
                        "t0_reco_no_isw", "t0_reio_no_isw", "t0_isw",
                        "t2_reco", "t2_reio"
                    ])

                # TODO SG: Maybe add further source functions if needed ...
                # if self.pt.has_source_delta_g:
                #     index_type_list.append(self.pt.index_tp_delta_g)
                #     names.append('delta_g')
                # if self.pt.has_source_theta_m:
                #     index_type_list.append(self.pt.index_tp_theta_m)
                #     names.append('theta_m')
                # if self.pt.has_source_phi:
                #     index_type_list.append(self.pt.index_tp_phi)
                #     names.append('phi')
                # if self.pt.has_source_phi_prime:
                #     index_type_list.append(self.pt.index_tp_phi_prime)
                #     names.append('phi_prime')
                # if self.pt.has_source_psi:
                #     index_type_list.append(self.pt.index_tp_psi)
                #     names.append('psi')
                timer.end("neural network initialization")  # 4e-5 sec

                timer.start("allocate numpy array of predictions")
                
                #try cython array init
                self.NN_prediction = cvarray(shape=(len(source_names),tau_size*k_NN_size),itemsize=sizeof(double),format='d')
                timer.end("allocate numpy array of predictions") # 1e-4 sec. Is at the 0.2% of runtime. But most likeli some space for improment ...

                # da is nen problem irgjendwoh. I have to init them as 0, otherwise something goes wrong somewhere
                for i in range(len(source_names)):
                    self.NN_prediction[i][:] = 0.0

                timer.start("get all sources")
                # NN prediction of the source functions. They are stored in the self.NN_prediction array.
                self.predictor.predict_many(source_names, np.asarray(tau_CLASS), self.NN_prediction)
                timer.end("get all sources") #0.07 sec

                timer.start("overwrite k array")
                # Copy k values from NN
                free(self.pt.k[index_md])
                self.pt.k[index_md] = <double*>malloc(k_NN_size * sizeof(double))
                for i_k in range(k_NN_size):
                    self.pt.k[index_md][i_k] = k_NN[i_k]
   
                self.pt.k_min = k_NN[0]
                self.pt.k_max = k_NN[-1]
                self.pt.k_size[index_md] = k_NN_size
                # the following part determines k_max_cl similar to the c function perturb_get_k_list in the source/perturbations.c file of ClassFull
                #first value
                if self.ba.sgnK == 0:
                    #K<0 (flat)  : start close to zero */
                    k_min=self.pr.k_min_tau0/self.ba.conformal_age
                elif self.ba.sgnK == -1:
                    # K<0 (open)  : start close to sqrt(-K)
                    # (in transfer modules, for scalars, this will correspond to q close to zero
                    # for vectors and tensors, this value is even smaller than the minimum necessary value) */
                    k_min=np.sqrt(-self.ba.K+pow(self.pr.k_min_tau0/self.ba.conformal_age/self.th.angular_rescaling,2))
                elif self.ba.sgnK == 1:
                    # K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
                    k_min = np.sqrt((8.-1.e-4)*self.ba.K)
                #- --> find k_max (as well as k_max_cmb[ppt->index_md_scalars], k_max_cl[ppt->index_md_scalars]) */
                k_max_cmb = k_min
                k_max_cl = k_min
                if self.pt.has_cls == _TRUE_:
                    # find k_max_cmb[ppt->index_md_scalars] : */
                    # choose a k_max_cmb[ppt->index_md_scalars] corresponding to a wavelength on the last
                    # scattering surface seen today under an angle smaller than
                    # pi/lmax: this is equivalent to
                    # k_max_cl[ppt->index_md_scalars]*[comvoving.ang.diameter.distance] > l_max */
                    k_max_cmb = self.pr.k_max_tau0_over_l_max*self.pt.l_scalar_max/self.ba.conformal_age/self.th.angular_rescaling
                    k_max_cl = k_max_cmb
                # find k_max_cl[ppt->index_md_scalars] : */
                # if we need density/lensing Cl's, we must impose a stronger condition,
                # such that the minimum wavelength on the shell corresponding
                # to the center of smallest redshift bin is seen under an
                # angle smaller than pi/lmax. So we must multiply our previous
                # k_max_cl[ppt->index_md_scalars] by the ratio tau0/(tau0-tau[center of smallest
                # redshift bin]). Note that we could do the same with the
                # lensing potential if we needed a very precise C_l^phi-phi at
                if self.pt.has_cl_number_count == _TRUE_ or self.pt.has_cl_lensing_potential == _TRUE_:
                    if background_tau_of_z(&self.ba,self.pt.selection_mean[0],&tau1) == _FAILURE_:
                        raise CosmoSevereError(self.ba.error_message)
                    k_max_cl = max(k_max_cl,self.pr.k_max_tau0_over_l_max*self.pt.l_lss_max/(self.ba.conformal_age-tau1))
                    # to be very accurate we should use angular diameter distance to given redshift instead of comoving radius: would implement corrections depending on curvature
                
                # end of c part 
                
                k_max_cl_idx = np.searchsorted(k_NN, k_max_cl)
                # self.pt.k_size_cl[index_md] = k_NN_size
                self.pt.k_size_cl[index_md] = k_max_cl_idx

                _k_max_dbg = self.pt.k[index_md][self.pt.k_size_cl[index_md] - 1]
                timer.end("overwrite k array") # 5e-5 sec

                # copy the predictor timestamps to the classy one
                for key, value in self.predictor.times.items():
                    timer[key] = value

                # copy the predictor timestamps to the classy one     #takes about 1e-5 sec
                for key, value in self.predictor.time_prediction_per_network.items():
                    timer["indiv. network: '{}'".format(key)] = value

                #allocate all source functions
                timer.start("allocate unused source functions")
                #define source size
                self.NN_source_size = k_NN_size * tau_size * sizeof(double)
                for index_type in range(tp_size):
                    # We need to allocate these, such that they can be dealocated by the perturb_free function later on
                    self.pt.sources[index_md][index_ic * tp_size + index_type] = <double*> malloc(sizeof(double))
                timer.end("allocate unused source functions") #takes about 1e-6 sec

                timer.start("overwrite source functions")
                start = time.time()
                for i, index_tp_x in enumerate(NN_requested_index_types):
                    self.overwrite_source_function(
                            index_md, index_ic,
                            index_tp_x,
                            k_NN_size, tau_size, 
                            self.NN_prediction[i, :]
                            )
                timer.end("overwrite source functions") # this takes 3e-6 sec
                timer.end("neural network complete")

            timer.end("perturb")

            if post_perturb_callback:
                post_perturb_callback(self)


        if "primordial" in level:
            timer.start("primordial")
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            self.ncp.add("primordial")
            timer.end("primordial")


        if "fourier" in level:
            timer.start("fourier")
            if fourier_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.fo) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.fo.error_message)
            self.ncp.add("fourier")
            timer.end("fourier")

        if "transfer" in level:
            timer.start("transfer")
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.fo), &(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)
            self.ncp.add("transfer")
            timer.end("transfer")

        if "harmonic" in level:
            timer.start("harmonic")
            if harmonic_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.fo), &(self.tr),
                            &(self.hr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.hr.error_message)
            self.ncp.add("harmonic")
            timer.end("harmonic")

        if "lensing" in level:
            timer.start("lensing")
            if lensing_init(&(self.pr), &(self.pt), &(self.hr),
                            &(self.fo), &(self.le)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.le.error_message)
            self.ncp.add("lensing")
            timer.end("lensing")


        if "distortions" in level:
            timer.start("distortions")
            if distortions_init(&(self.pr), &(self.ba), &(self.th),
                                &(self.pt), &(self.pm), &(self.sd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sd.error_message)
            self.ncp.add("distortions")
            timer.end("distortions")

        self.computed = True

        timer.end("compute")

        if performance_report is not None:
            performance_report.update(timer.times)

        # At this point, the cosmological instance contains everything needed. The
        # following functions are only to output the desired numbers
        return 

    def raw_cl(self, lmax=-1, nofail=False):
        """
        raw_cl(lmax=-1, nofail=False)

        Return a dictionary of the primary C_l

        Parameters
        ----------
        lmax : int, optional
                Define the maximum l for which the C_l will be returned
                (inclusively). This number will be checked against the maximum l
                at which they were actually computed by CLASS, and an error will
                be raised if the desired lmax is bigger than what CLASS can
                give.
        nofail: bool, optional
                Check and enforce the computation of the harmonic module
                beforehand, with the desired lmax.

        Returns
        -------
        cl : dict
                Dictionary that contains the power spectrum for 'tt', 'te', etc... The
                index associated with each is defined wrt. Class convention, and are non
                important from the python point of view. It also returns now the
                ell array.
        """
        cdef int lmaxR
        cdef double *rcl = <double*> calloc(self.hr.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md[index_md] = <double*> calloc(self.hr.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.hr.ct_size*self.hr.ic_ic_size[index_md], sizeof(double))

        # Define a list of integers, refering to the flags and indices of each
        # possible output Cl. It allows for a clear and concise way of looping
        # over them, checking if they are defined or not.
        has_flags = [
            (self.hr.has_tt, self.hr.index_ct_tt, 'tt'),
            (self.hr.has_ee, self.hr.index_ct_ee, 'ee'),
            (self.hr.has_te, self.hr.index_ct_te, 'te'),
            (self.hr.has_bb, self.hr.index_ct_bb, 'bb'),
            (self.hr.has_pp, self.hr.index_ct_pp, 'pp'),
            (self.hr.has_tp, self.hr.index_ct_tp, 'tp'),]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)

        if not spectra:
            raise CosmoSevereError("No Cl computed")
        lmaxR = self.hr.l_max_tot
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        # Initialise all the needed Cls arrays
        cl = {}
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)

        # Recover for each ell the information from CLASS
        for ell from 2<=ell<lmax+1:
            if harmonic_cl_at_l(&self.hr, ell, rcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.hr.error_message)
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = rcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(rcl)
        for index_md in range(self.hr.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def lensed_cl(self, lmax=-1,nofail=False):
        """
        lensed_cl(lmax=-1, nofail=False)

        Return a dictionary of the lensed C_l, computed by CLASS, without the
        density C_ls. They must be asked separately with the function aptly
        named density_cl

        Parameters
        ----------
        lmax : int, optional
                Define the maximum l for which the C_l will be returned (inclusively)
        nofail: bool, optional
                Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : dict
                Dictionary that contains the power spectrum for 'tt', 'te', etc... The
                index associated with each is defined wrt. Class convention, and are non
                important from the python point of view.
        """
        cdef int lmaxR
        cdef double *lcl = <double*> calloc(self.le.lt_size,sizeof(double))

        # Define a list of integers, refering to the flags and indices of each
        # possible output Cl. It allows for a clear and concise way of looping
        # over them, checking if they are defined or not.
        has_flags = [
            (self.le.has_tt, self.le.index_lt_tt, 'tt'),
            (self.le.has_ee, self.le.index_lt_ee, 'ee'),
            (self.le.has_te, self.le.index_lt_te, 'te'),
            (self.le.has_bb, self.le.index_lt_bb, 'bb'),
            (self.le.has_pp, self.le.index_lt_pp, 'pp'),
            (self.le.has_tp, self.le.index_lt_tp, 'tp'),]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)

        if not spectra:
            raise CosmoSevereError("No lensed Cl computed")
        lmaxR = self.le.l_lensed_max

        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}
        # Simple Cls, for temperature and polarisation, are not so big in size
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)
        for ell from 2<=ell<lmax+1:
            if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
                raise CosmoSevereError(self.le.error_message)
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = lcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(lcl)
        return cl

    def density_cl(self, lmax=-1, nofail=False):
        """
        density_cl(lmax=-1, nofail=False)

        Return a dictionary of the primary C_l for the matter

        Parameters
        ----------
        lmax : int, optional
            Define the maximum l for which the C_l will be returned (inclusively)
        nofail: bool, optional
            Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : numpy array of numpy.ndarrays
            Array that contains the list (in this order) of self correlation of
            1st bin, then successive correlations (set by non_diagonal) to the
            following bins, then self correlation of 2nd bin, etc. The array
            starts at index_ct_dd.
        """
        cdef int lmaxR
        cdef double *dcl = <double*> calloc(self.hr.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md[index_md] = <double*> calloc(self.hr.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.hr.ct_size*self.hr.ic_ic_size[index_md], sizeof(double))

        lmaxR = self.pt.l_lss_max
        has_flags = [
            (self.hr.has_dd, self.hr.index_ct_dd, 'dd'),
            (self.hr.has_td, self.hr.index_ct_td, 'td'),
            (self.hr.has_ll, self.hr.index_ct_ll, 'll'),
            (self.hr.has_dl, self.hr.index_ct_dl, 'dl'),
            (self.hr.has_tl, self.hr.index_ct_tl, 'tl')]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)
                l_max_flag = self.hr.l_max_ct[self.hr.index_md_scalars][index]
                if l_max_flag < lmax and lmax > 0:
                    raise CosmoSevereError(
                        "the %s spectrum was computed until l=%i " % (
                            name.upper(), l_max_flag) +
                        "but you asked a l=%i" % lmax)

        if not spectra:
            raise CosmoSevereError("No density Cl computed")
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_lss",lmax)
                self._pars_check("output",'nCl')
                self.compute()
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}

        # For density Cls, the size is bigger (different redshfit bins)
        # computes the size, given the number of correlations needed to be computed
        size = int((self.hr.d_size*(self.hr.d_size+1)-(self.hr.d_size-self.hr.non_diag)*
                (self.hr.d_size-1-self.hr.non_diag))/2);
        for elem in ['dd', 'll', 'dl']:
            if elem in spectra:
                cl[elem] = {}
                for index in range(size):
                    cl[elem][index] = np.zeros(
                        lmax+1, dtype=np.double)
        for elem in ['td', 'tl']:
            if elem in spectra:
                cl[elem] = np.zeros(lmax+1, dtype=np.double)

        for ell from 2<=ell<lmax+1:
            if harmonic_cl_at_l(&self.hr, ell, dcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.hr.error_message)
            if 'dd' in spectra:
                for index in range(size):
                    cl['dd'][index][ell] = dcl[self.hr.index_ct_dd+index]
            if 'll' in spectra:
                for index in range(size):
                    cl['ll'][index][ell] = dcl[self.hr.index_ct_ll+index]
            if 'dl' in spectra:
                for index in range(size):
                    cl['dl'][index][ell] = dcl[self.hr.index_ct_dl+index]
            if 'td' in spectra:
                cl['td'][ell] = dcl[self.hr.index_ct_td]
            if 'tl' in spectra:
                cl['tl'][ell] = dcl[self.hr.index_ct_tl]
        cl['ell'] = np.arange(lmax+1)

        free(dcl)
        for index_md in range(self.hr.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def z_of_r (self,z_array):
        cdef int last_index=0 #junk
        cdef double * pvecback
        r = np.zeros(len(z_array),'float64')
        dzdr = np.zeros(len(z_array),'float64')

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        i = 0
        for redshift in z_array:

            if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # store r
            r[i] = pvecback[self.ba.index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[self.ba.index_bg_H]

            i += 1

        free(pvecback)
        return r[:],dzdr[:]

    def luminosity_distance(self, z):
        """
        luminosity_distance(z)
        """
        cdef int last_index = 0  # junk
        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba, z, long_info,
                inter_normal, &last_index, pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        lum_distance = pvecback[self.ba.index_bg_lum_distance]
        free(pvecback)
        return lum_distance

    # Gives the total matter pk for a given (k,z)
    def pk(self,double k,double z):
        """
        Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if (self.fo.method == nl_none):
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_m,&pk,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)
        else:
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_nonlinear,k,z,self.fo.index_pk_m,&pk,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return pk

    # Gives the total matter pk for a given (z)
    def pk_at_z(self,double z):
        """

        Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """

        k = np.zeros(self.fo.k_size,'float64')
        pk = np.zeros(self.fo.k_size,'float64')

        pk_out = <double*> malloc(self.fo.k_size*sizeof(double))
        ln_pk_ic = <double*> malloc(self.fo.ic_ic_size*self.fo.k_size*sizeof(double))

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")
        # logarithmic = 1
        if (self.fo.method == nl_none):
            if fourier_pk_at_z(&self.ba,&self.fo,linear,pk_linear,z,self.fo.index_pk_m,pk_out,ln_pk_ic)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)
        else:
            if fourier_pk_at_z(&self.ba,&self.fo,linear,pk_nonlinear,z,self.fo.index_pk_m,pk_out,ln_pk_ic)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        for i in range(self.fo.k_size):
            pk[i] = pk_out[i]
            k[i] = np.exp(self.fo.ln_k[i])

        free(ln_pk_ic)
        free(pk_out)

        return pk, k


    # Gives the total matter pk for a given (z)
    def pk_cb_at_z(self,double z):
        """

        Gives the cdm + b pk (in Mpc**3) for a given z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """

        k = np.zeros(self.fo.k_size,'float64')
        pk_cb = np.zeros(self.fo.k_size,'float64')

        pk_cb_out = <double*> malloc(self.fo.k_size*sizeof(double))
        ln_pk_ic = <double*> malloc(self.fo.ic_ic_size*self.fo.k_size*sizeof(double))

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")
        # logarithmic = 1
        if (self.fo.method == nl_none):
            if fourier_pk_at_z(&self.ba,&self.fo,linear,pk_linear,z,self.fo.index_pk_cb,pk_cb_out,ln_pk_ic)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)
        else:
            if fourier_pk_at_z(&self.ba,&self.fo,linear,pk_nonlinear,z,self.fo.index_pk_cb,pk_cb_out,ln_pk_ic)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        for i in range(self.fo.k_size):
            pk_cb[i] = pk_cb_out[i]
            k[i] = np.exp(self.fo.ln_k[i])

        free(ln_pk_ic)
        free(pk_cb_out)

        return pk_cb, k

    # Gives the cdm+b pk for a given (k,z)
    def pk_cb(self,double k,double z):
        """
        Gives the cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_cb

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")
        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        if (self.fo.method == nl_none):
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)
        else:
            if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_nonlinear,k,z,self.fo.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return pk_cb

    # Gives the total matter pk for a given (k,z)
    def pk_lin(self,double k,double z):
        """
        Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_m,&pk_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_lin

    # Gives the cdm+b pk for a given (k,z)
    def pk_cb_lin(self,double k,double z):
        """
        Gives the linear cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_cb_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_cb,&pk_cb_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_cb_lin

    def get_pk(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    def get_pk_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    # [NS] :: TODO :: check optimization
    def get_pk_all(self, k, z, nonlinear = True, cdmbar = False, z_axis_in_k_arr = 0):
        """ General function to get the P(k,z) for ARBITRARY shapes of k,z
            Additionally, it includes the functionality of selecting wether to use the non-linear parts or not,
            and wether to use the cdm baryon power spectrum only
            For Multi-Dimensional k-arrays, it assumes that one of the dimensions is the z-axis
            This is handled by the z_axis_in_k_arr integer, as described in the source code """
        # z_axis_in_k_arr specifies the integer position of the z_axis wihtin the n-dimensional k_arr
        # Example: 1-d k_array -> z_axis_in_k_arr = 0
        # Example: 3-d k_array with z_axis being the first axis -> z_axis_in_k_arr = 0
        # Example: 3-d k_array with z_axis being the last axis  -> z_axis_in_k_arr = 2


        # 1) Select the correct function
        if nonlinear:
            if cdmbar and not (self.ba.Omega0_ncdm_tot == 0.):
                pk_function = self.pk_cb
            else:
                pk_function = self.pk
        else:
            if cdmbar and not (self.ba.Omega0_ncdm_tot == 0.):
                pk_function = self.pk_cb_lin
            else:
                pk_function = self.pk_lin

        # 2) Check if z array, or z value
        if not isinstance(z,(list,np.ndarray)):
            # Only single z value was passed -> k could still be an array of arbitrary dimension
            if not isinstance(k,(list,np.ndarray)):
                # Only single z value AND only single k value -> just return a value
                # This iterates over ALL remaining dimensions
                return pk_function(k,z)
            else:
                k_arr = np.array(k)
                out_pk = np.empty(np.shape(k_arr))
                iterator = np.nditer(k_arr,flags=['multi_index'])
                while not iterator.finished:
                    out_pk[iterator.multi_index] = pk_function(iterator[0],z)
                    iterator.iternext()
                # This iterates over ALL remaining dimensions
                #for index_k in range(k_arr.shape[-1]):
                #    out_pk[...,index_k] = pk_function(k_arr[...,index_k],z)
                return out_pk

        # 3) An array of z values was passed
        k_arr = np.array(k)
        z_arr = np.array(z)
        if( z_arr.ndim != 1 ):
            raise CosmoSevereError("Can only parse one-dimensional z-arrays, not multi-dimensional")

        if( k_arr.ndim > 1 ):
            # 3.1) If there is a multi-dimensional k-array of EQUAL lenghts
            out_pk = np.empty(np.shape(k_arr))
            # Bring the z_axis to the front
            k_arr = np.moveaxis(k_arr, z_axis_in_k_arr, 0)
            out_pk = np.moveaxis(out_pk, z_axis_in_k_arr, 0)
            if( len(k_arr) != len(z_arr) ):
                raise CosmoSevereError("Mismatching array lengths of the z-array")
            for index_z in range(len(z_arr)):
                iterator = np.nditer(k_arr[index_z],flags=['multi_index'])
                while not iterator.finished:
                    out_pk[index_z][iterator.multi_index] = pk_function(iterator[0],z[index_z])
                    iterator.iternext()
                # This iterates over ALL remaining dimensions
                #for index_k in range(k_arr[index_z].shape[-1]):
                #    out_pk[index_z][...,index_k] = pk_function(k_arr[index_z][...,index_k],z_arr[index_z])
            # Move the z_axis back into position
            k_arr = np.moveaxis(k_arr, 0, z_axis_in_k_arr)
            out_pk = np.moveaxis(out_pk, 0, z_axis_in_k_arr)
            return out_pk
        else:
            # 3.2) If there is a multi-dimensional k-array of UNEQUAL lenghts
            if isinstance(k_arr[0],(list,np.ndarray)):
                # A very special thing happened: The user passed a k array with UNEQUAL lengths of k arrays for each z
                out_pk = []
                for index_z in range(len(z_arr)):
                    k_arr_at_z = np.array(k_arr[index_z])
                    out_pk_at_z = np.empty(np.shape(k_arr_at_z))
                    iterator = np.nditer(k_arr_at_z,flags=['multi_index'])
                    while not iterator.finished:
                        out_pk_at_z[iterator.multi_index] = pk_function(iterator[0],z[index_z])
                        iterator.iternext()
                    #for index_k in range(k_arr_at_z.shape[-1]):
                    #   out_pk_at_z[...,index_k] = pk_function(k_arr_at_z[...,index_k],z_arr[index_z])
                    out_pk.append(out_pk_at_z)
                return out_pk

            # 3.3) If there is a single-dimensional k-array
            # The user passed a z-array, but only a 1-d k array
            # Assume thus, that the k array should be reproduced for all z
            out_pk = np.empty((len(z_arr),len(k_arr)))
            for index_z in range(len(z_arr)):
                for index_k in range(len(k_arr)):
                    out_pk[index_z,index_k] = pk_function(k_arr[index_k],z_arr[index_z])
            return out_pk

    #################################
    # Gives a grid of values of matter and/or cb power spectrum, together with the vectors of corresponding k and z values
    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species = False, h_units=False):
        """
        Returns a grid of matter power spectrum values and the z and k
        at which it has been fully computed. Useful for creating interpolators.

        Parameters
        ----------
        nonlinear : bool
                Whether the returned power spectrum values are linear or non-linear (default)
        only_clustering_species : bool
                Whether the returned power spectrum is for galaxy clustering and excludes massive neutrinos, or always includes everything (default)
        h_units : bool
                Whether the units of k in output are h/Mpc or 1/Mpc (default)

        Returns
        -------
        pk : grid of power spectrum values, pk[index_k,index_z]
        k : vector of k values, k[index_k] (in units of 1/Mpc by default, or h/Mpc when setting h_units to True)
        z : vector of z values, z[index_z]
        """

        cdef np.ndarray[DTYPE_t,ndim=2] pk = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.fo.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.fo.k_size_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.fo.ln_tau_size-self.fo.index_ln_tau_pk),'float64')
        cdef int index_k, index_tau, index_pk
        cdef double z_max_nonlinear, z_max_requested

        # consistency checks

        if self.fo.has_pk_matter == False:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations; add 'mPk' in 'output'")

        if nonlinear == True and self.fo.method == nl_none:
            raise CosmoSevereError("You ask classy to return an array of nonlinear P(k,z) values, but the input parameters sent to CLASS did not require any non-linear P(k,z) calculations; add e.g. 'halofit' or 'HMcode' in 'nonlinear'")

        # check wich type of P(k) to return (total or clustering only, i.e. without massive neutrino contribution)
        if (only_clustering_species == True):
            index_pk = self.fo.index_pk_cluster
        else:
            index_pk = self.fo.index_pk_total

        # get list of redshifts
        # the ln(times) of interest are stored in self.fo.ln_tau[index_tau]
        # with index_tau in the range:
        # self.fo.index_ln_tau_pk <= index_tau < self.fo.ln_tau_size-self.fo.index_ln_tau_pk
        #
        # If the user forgot to ask for z_max_pk>0, then self.fo.ln_tau_size == 1 and self.fo.index_ln_tau_pk==0

        if self.fo.ln_tau_size == 1:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in xrange(self.fo.ln_tau_size-self.fo.index_ln_tau_pk):
                if index_tau == self.fo.ln_tau_size-self.fo.index_ln_tau_pk-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.fo.ln_tau[index_tau+self.fo.index_ln_tau_pk]))

        # check consitency of the list of redshifts

        if nonlinear == True:
            # Check highest value of z at which nl corrections could be computed.
            # In the table tau_sampling it corresponds to index: self.fo.index_tau_min_n
            z_max_nonlinear = self.z_of_tau(self.fo.tau[self.fo.index_tau_min_nl])

            # Check highest value of z in the requested output.
            # In the table tau_sampling it corresponds to index: self.fo.tau_size - self.fo.ln_tau_size + self.fo.index_ln_tau_pk
            z_max_requested = z[0]

            # The first z must be larger or equal to the second one, that is,
            # the first index must be smaller or equal to the second one.
            # If not, raise and error.

            if (self.fo.index_tau_min_nl > self.fo.tau_size - self.fo.ln_tau_size + self.fo.index_ln_tau_pk):
                raise CosmoSevereError("get_pk_and_k_and_z() is trying to return P(k,z) up to z_max=%e (to encompass your requested maximum value of z); but the input parameters sent to CLASS (in particular ppr->nonlinear_min_k_max=%e) were such that the non-linear P(k,z) could only be consistently computed up to z=%e; increase the precision parameter 'nonlinear_min_k_max', or decrease your requested z_max"%(z_max_requested,self.pr.nonlinear_min_k_max,z_max_nonlinear))

        # get list of k

        if h_units:
            units=1./self.ba.h
        else:
            units=1

        for index_k in xrange(self.fo.k_size_pk):
            k[index_k] = self.fo.k[index_k]*units

        # get P(k,z) array

        for index_tau in xrange(self.fo.ln_tau_size-self.fo.index_ln_tau_pk):
            for index_k in xrange(self.fo.k_size_pk):
                if nonlinear == True:
                    pk[index_k, index_tau] = np.exp(self.fo.ln_pk_nl[index_pk][(index_tau+self.fo.index_ln_tau_pk) * self.fo.k_size + index_k])
                else:
                    pk[index_k, index_tau] = np.exp(self.fo.ln_pk_l[index_pk][(index_tau+self.fo.index_ln_tau_pk) * self.fo.k_size + index_k])

        return pk, k, z

    #################################
    # Gives a grid of each transfer functions arranged in a dictionary, together with the vectors of corresponding k and z values
    def get_transfer_and_k_and_z(self, output_format='class', h_units=False):
        """
        Returns a dictionary of grids of density and/or velocity transfer function values and the z and k at which it has been fully computed.
        Useful for creating interpolators.
        When setting CLASS input parameters, include at least one of 'dTk' (for density transfer functions) or 'vTk' (for velocity transfer functions).
        Following the default output_format='class', all transfer functions will be normalised to 'curvature R=1' at initial time
        (and not 'curvature R = -1/k^2' like in CAMB).
        You may switch to output_format='camb' for the CAMB definition and normalisation of transfer functions.
        (Then, 'dTk' must be in the input: the CAMB format only outputs density transfer functions).
        When sticking to output_format='class', you also get the newtonian metric fluctuations phi and psi.
        If you set the CLASS input parameter 'extra_metric_transfer_functions' to 'yes',
        you get additional metric fluctuations in the synchronous and N-body gauges.

        Parameters
        ----------
        output_format  : ('class' or 'camb')
                Format transfer functions according to CLASS (default) or CAMB
        h_units : bool
                Whether the units of k in output are h/Mpc or 1/Mpc (default)

        Returns
        -------
        tk : dictionary containing all transfer functions.
                For instance, the grid of values of 'd_c' (= delta_cdm) is available in tk['d_c']
                All these grids have indices [index_k,index,z], for instance tk['d_c'][index_k,index,z]
        k : vector of k values (in units of 1/Mpc by default, or h/Mpc when setting h_units to True)
        z : vector of z values
        """
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.pt.k_size_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.pt.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef int index_k, index_tau
        cdef char * titles
        cdef double * data
        cdef file_format outf

        # consistency checks
        if (self.pt.has_density_transfers == False) and (self.pt.has_velocity_transfers == False):
            raise CosmoSevereError("You ask classy to return transfer functions, but the input parameters sent to CLASS did not require any T(k,z) calculations; add 'dTk' and/or 'vTk' in 'output'")

        index_md = self.pt.index_md_scalars;

        if (self.pt.ic_size[index_md] > 1):
            raise CosmoSevereError("For simplicity, get_transfer_and_k_and_z() has been written assuming only adiabatic initial conditions. You need to write the generalisation to cases with multiple initial conditions.")

        # check out put format
        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        # check name and number of trnasfer functions computed ghy CLASS

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if perturbations_output_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)

        # get list of redshifts
        # the ln(times) of interest are stored in self.fo.ln_tau[index_tau]
        # with index_tau in the range:
        # self.fo.index_ln_tau_pk <= index_tau < self.fo.ln_tau_size-self.fo.index_ln_tau_pk
        #
        # If the user forgot to ask for z_max_pk>0, then self.fo.ln_tau_size == 1 and self.fo.index_ln_tau_pk==0

        z_size = self.pt.ln_tau_size-self.pt.index_ln_tau_pk
        if z_size == 1:
            raise CosmoSevereError("You ask classy to return an array of T_x(k,z) values, but the input parameters sent to CLASS did not require any transfer function calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in xrange(z_size):
                if index_tau == z_size-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.pt.ln_tau[index_tau+self.pt.index_ln_tau_pk]))

        # get list of k

        if h_units:
            units=1./self.ba.h
        else:
            units=1

        k_size = self.pt.k_size_pk
        for index_k in xrange(k_size):
            k[index_k] = self.pt.k[index_md][index_k]*units

        # create output dictionary

        tk = {}
        for index_type,name in enumerate(names):
            if index_type > 0:
                tk[name] = np.zeros((k_size, z_size),'float64')

        # allocate the vector in wich the transfer functions will be stored temporarily for all k and types at a given z
        data = <double*>malloc(sizeof(double)*number_of_titles*self.pt.k_size[index_md])

        # get T(k,z) array

        for index_tau in xrange(z_size):

            if perturbations_output_data_at_index_tau(&self.ba, &self.pt, outf, index_tau+self.pt.index_ln_tau_pk, number_of_titles, data)==_FAILURE_:
                raise CosmoSevereError(self.pt.error_message)

            for index_type,name in enumerate(names):
                if index_type > 0:
                    for index_k in xrange(k_size):
                        tk[name][index_k, index_tau] = data[index_k*number_of_titles+index_type]

        free(data)
        return tk, k, z

    #################################
    # Gives a grid of values of the power spectrum of the quantity [k^2*(phi+psi)/2], where (phi+psi)/2 is the Weyl potential, together with the vectors of corresponding k and z values
    def get_Weyl_pk_and_k_and_z(self, nonlinear=False, h_units=False):
        """
        Returns a grid of Weyl potential (phi+psi) power spectrum values and the z and k
        at which it has been fully computed. Useful for creating interpolators.
        Note that this function just calls get_pk_and_k_and_z and corrects the output
        by the ratio of transfer functions [(phi+psi)/d_m]^2.

        Parameters
        ----------
        nonlinear : bool
                Whether the returned power spectrum values are linear or non-linear (default)
        h_units : bool
                Whether the units of k in output are h/Mpc or 1/Mpc (default)

        Returns
        -------
        Weyl_pk : grid of Weyl potential (phi+psi) spectrum values, Weyl_pk[index_k,index_z]
        k : vector of k values, k[index_k] (in units of 1/Mpc by default, or h/Mpc when setting h_units to True)
        z : vector of z values, z[index_z]
        """
        cdef np.ndarray[DTYPE_t,ndim=2] pk = np.zeros((self.fo.k_size_pk,self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] k4 = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] phi = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] psi = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] d_m = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] Weyl_pk = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size-self.pt.index_ln_tau_pk),'float64')

        cdef bint input_nonlinear = nonlinear
        cdef bint input_h_units = h_units

        cdef int index_z

        # get total matter power spectrum
        pk, k, z = self.get_pk_and_k_and_z(nonlinear=input_nonlinear, only_clustering_species = False, h_units=input_h_units)

        # get transfer functions
        tk_and_k_and_z = {}
        tk_and_k_and_z, k, z = self.get_transfer_and_k_and_z(output_format='class',h_units=input_h_units)
        phi = tk_and_k_and_z['phi']
        psi = tk_and_k_and_z['psi']
        d_m = tk_and_k_and_z['d_m']

        # get an array containing k**4 (same for all redshifts)
        for index_z in range(self.fo.ln_tau_size-self.pt.index_ln_tau_pk):
            k4[:,index_z] = k**4

        # rescale total matter power spectrum to get the Weyl power spectrum times k**4
        # (the latter factor is just a convention. Since there is a factor k**2 in the Poisson equation,
        # this rescaled Weyl spectrum has a shape similar to the matter power spectrum).
        Weyl_pk = pk * ((phi+psi)/2./d_m)**2 * k4

        return Weyl_pk, k, z

    #################################
    # Gives sigma(R,z) for a given (R,z)
    def sigma(self,double R,double z, h_units = False):
        """
        Gives sigma (total matter) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z).
         This is unless h_units is set to true, in which case R is the radius in units of Mpc/h,
         and R=8 corresponds to sigma8(z))

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma

        R_in_Mpc = (R if not h_units else R/self.ba.h)

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        if fourier_sigmas_at_z(&self.pr,&self.ba,&self.fo,R_in_Mpc,z,self.fo.index_pk_m,out_sigma,&sigma)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return sigma

    # Gives sigma_cb(R,z) for a given (R,z)
    def sigma_cb(self,double R,double z, h_units = False):
        """
        Gives sigma (cdm+b) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)
         This is unless h_units is set to true, in which case R is the radius in units of Mpc/h,
         and R=8 corresponds to sigma8(z))

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma_cb

        R_in_Mpc = (R if not h_units else R/self.ba.h)

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("sigma_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        # If necessary, convert R to units of Mpc
        if h_units:
            R /= self.ba.h

        if fourier_sigmas_at_z(&self.pr,&self.ba,&self.fo,R,z,self.fo.index_pk_cb,out_sigma,&sigma_cb)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return sigma_cb

    # Gives effective logarithmic slope of P_L(k,z) (total matter) for a given (k,z)
    def pk_tilt(self,double k,double z):
        """
        Gives effective logarithmic slope of P_L(k,z) (total matter) for a given k and z
        (k is the wavenumber in units of 1/Mpc, z is the redshift, the output is dimensionless)

        .. note::

            there is an additional check to verify whether output contains `mPk` and whether k is in the right range

        """
        cdef double pk_tilt

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get pk_tilt(k,z) you must add mPk to the list of outputs.")

        if (k < self.fo.k[1] or k > self.fo.k[self.fo.k_size-2]):
            raise CosmoSevereError("In order to get pk_tilt at k=%e 1/Mpc, you should compute P(k,z) in a wider range of k's"%k)

        if fourier_pk_tilt_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_total,&pk_tilt)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_tilt

    #calculates the hmcode window_function of the Navarrow Frenk White Profile
    def fourier_hmcode_window_nfw(self,double k,double rv,double c):
        """
        Gives window_nfw for a given wavevector k, virial radius rv and concentration c

        """
        cdef double window_nfw


        if fourier_hmcode_window_nfw(&self.fo,k,rv,c,&window_nfw)==_FAILURE_:
                 raise CosmoSevereError(self.hr.error_message)

        return window_nfw

    def age(self):
        self.compute(["background"])
        return self.ba.age

    def h(self):
        return self.ba.h

    def n_s(self):
        return self.pm.n_s

    def tau_reio(self):
        return self.th.tau_reio

    def Omega_m(self):
        return self.ba.Omega0_m

    def Omega_r(self):
        return self.ba.Omega0_r

    def theta_s_100(self):
        return 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)

    def theta_star_100(self):
        return 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)

    def omega_cdm(self):
        return self.ba.Omega0_cdm * self.ba.h * self.ba.h

    def Omega_Lambda(self):
        return self.ba.Omega0_lambda

    def Omega_g(self):
        return self.ba.Omega0_g

    def Omega_b(self):
        return self.ba.Omega0_b

    def omega_b(self):
        return self.ba.Omega0_b * self.ba.h * self.ba.h

    def Neff(self):
        return self.ba.Neff

    def k_eq(self):
        self.compute(["background"])
        return self.ba.a_eq*self.ba.H_eq

    def z_eq(self):
        self.compute(["background"])
        return 1./self.ba.a_eq-1.

    def sigma8(self):
        self.compute(["fourier"])
        return self.fo.sigma8[self.fo.index_pk_m]

    #def neff(self):
    #    self.compute(["harmonic"])
    #    return self.hr.neff

    def sigma8_cb(self):
        self.compute(["fourier"])
        return self.fo.sigma8[self.fo.index_pk_cb]

    def rs_drag(self):
        self.compute(["thermodynamics"])
        return self.th.rs_d

    def rs_drag_nn(self):
        """
        Same as `self.rs_drag()`, but doesn't invoke `self.compute()`.
        The reason is the following: If NNs are enabled, `self.compute()`
        will call the NN code during the perturbation module; the NN code will
        call `rs_drag()`, which in turn will call `self.compute(["thermodynamics"])`;
        however,  since `self.ready` is not yet set to `True` during the perturbation
        module, this will recompute the thermodynamics (and waste time).
        For this reason, this method assumes that thermodynamics has been run already
        WITHOUT checking the `self.ready` flag.
        """
        return self.th.rs_d

    def z_reio(self):
        self.compute(["thermodynamics"])
        return self.th.z_reio

    def angular_distance(self, z):
        """
        angular_distance(z)

        Return the angular diameter distance (exactly, the quantity defined by Class
        as index_bg_ang_distance in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return D_A

    #################################
    # Get angular diameter distance of object at z2 as seen by observer at z1,
    def angular_distance_from_to(self, z1, z2):
        """
        angular_distance_from_to(z)

        Return the angular diameter distance of object at z2 as seen by observer at z1,
        that is, sin_K((chi2-chi1)*np.sqrt(|k|))/np.sqrt(|k|)/(1+z2).
        If z1>z2 returns zero.

        Parameters
        ----------
        z1 : float
                Observer redshift
        z2 : float
                Source redshift

        Returns
        -------
        d_A(z1,z2) in Mpc
        """
        cdef int last_index #junk
        cdef double * pvecback

        if z1>=z2:
            return 0.

        else:
            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z1,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # This is the comoving distance to object at z1
            chi1 = pvecback[self.ba.index_bg_conf_distance]

            if background_at_z(&self.ba,z2,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # This is the comoving distance to object at z2
            chi2 = pvecback[self.ba.index_bg_conf_distance]

            free(pvecback)

            if self.ba.K == 0:
                return (chi2-chi1)/(1+z2)
            elif self.ba.K > 0:
                return np.sin(np.sqrt(self.ba.K)*(chi2-chi1))/np.sqrt(self.ba.K)/(1+z2)
            elif self.ba.K < 0:
                return np.sinh(np.sqrt(-self.ba.K)*(chi2-chi1))/np.sqrt(-self.ba.K)/(1+z2)

    def comoving_distance(self, z):
        """
        comoving_distance(z)

        Return the comoving distance

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        r = pvecback[self.ba.index_bg_conf_distance]

        free(pvecback)

        return r

    def scale_independent_growth_factor(self, z):
        """
        scale_independent_growth_factor(z)
        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)
        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D = pvecback[self.ba.index_bg_D]

        free(pvecback)

        return D
        
    def scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale independent growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        f = pvecback[self.ba.index_bg_f]

        free(pvecback)

        return f

    #################################
    def scale_dependent_growth_factor_f(self, k, z, h_units=False, nonlinear=False, Nz=20):
        """
        scale_dependent_growth_factor_f(k,z)

        Return the scale dependent growth factor
        f(z)= 1/2 * [d ln P(k,a) / d ln a]
            = - 0.5 * (1+z) * [d ln P(k,z) / d z]
        where P(k,z) is the total matter power spectrum

        Parameters
        ----------
        z : float
                Desired redshift
        k : float
                Desired wavenumber in 1/Mpc (if h_units=False) or h/Mpc (if h_units=True)
        """

        # build array of z values at wich P(k,z) was pre-computed by class (for numerical derivative)
        # check that P(k,z) was stored at different zs
        if self.fo.ln_tau_size > 1:
            # check that input z is in stored range
            z_max = self.z_of_tau(np.exp(self.fo.ln_tau[0]))
            if (z<0) or (z>z_max):
                raise CosmoSevereError("You asked for f(k,z) at a redshift %e outside of the computed range [0,%e]"%(z,z_max))
            # create array of zs in growing z order (decreasing tau order)
            z_array = np.empty(self.fo.ln_tau_size)
            for i in xrange(self.fo.ln_tau_size):
                z_array[i] = self.z_of_tau(np.exp(self.fo.ln_tau[self.fo.ln_tau_size-1-i]))
            # due to interpolation errors, z_array[0] might be different from zero; round it to its true value.
            z_array[0]=0.
        else:
            raise CosmoSevereError("You asked for the scale-dependent growth factor: this requires numerical derivation of P(k,z) w.r.t z, and thus passing a non-zero input parameter z_max_pk")

        # if needed, convert k to units of 1/Mpc
        if h_units:
            k = k*self.ba.h

        # Allocate an array of P(k,z[...]) values
        Pk_array = np.empty_like(z_array)

        # Choose whether to use .pk() or .pk_lin()
        # The linear pk is in .pk_lin if nonlinear corrections have been computed, in .pk otherwise
        # The non-linear pk is in .pk if nonlinear corrections have been computed
        if nonlinear == False:
            if self.fo.method == nl_none:
                use_pk_lin = False
            else:
                use_pk_lin = True
        else:
            if self.fo.method == nl_none:
                raise CosmoSevereError("You asked for the scale-dependent growth factor of non-linear matter fluctuations, but you did not ask for non-linear calculations at all")
            else:
                use_pk_lin = False

        # Get P(k,z) and array P(k,z[...])
        if use_pk_lin == False:
            Pk = self.pk(k,z)
            for iz, zval in enumerate(z_array):
                Pk_array[iz] = self.pk(k,zval)
        else:
            Pk = self.pk_lin(k,z)
            for iz, zval in enumerate(z_array):
                Pk_array[iz] = self.pk_lin(k,zval)

        # Compute derivative (d ln P / d ln z)
        dPkdz = UnivariateSpline(z_array,Pk_array,s=0).derivative()(z)

        # Compute growth factor f
        f = -0.5*(1+z)*dPkdz/Pk

        return f

    #################################
    # gives f(z)*sigma8(z) where f(z) is the scale-independent growth factor
    def scale_independent_f_sigma8(self, z):
        """
        scale_independent_f_sigma8(z)

        Return the scale independent growth factor f(z) multiplied by sigma8(z)

        Parameters
        ----------
        z : float
                Desired redshift

        Returns
        -------
        f(z)*sigma8(z) (dimensionless)
        """
        return self.scale_independent_growth_factor_f(z)*self.sigma(8,z,h_units=True)

    #################################
    # gives an estimation of f(z)*sigma8(z) at the scale of 8 h/Mpc, computed as (d sigma8/d ln a)
    def effective_f_sigma8(self, z, z_step=0.1):
        """
        effective_f_sigma8(z)

        Returns the time derivative of sigma8(z) computed as (d sigma8/d ln a)

        Parameters
        ----------
        z : float
                Desired redshift
        z_step : float
                Default step used for the numerical two-sided derivative. For z < z_step the step is reduced progressively down to z_step/10 while sticking to a double-sided derivative. For z< z_step/10 a single-sided derivative is used instead.

        Returns
        -------
        (d ln sigma8/d ln a)(z) (dimensionless)
        """

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)

        # if possible, use two-sided derivative with default value of z_step
        if z >= z_step:
            return (self.sigma(8,z-z_step,h_units=True)-self.sigma(8,z+z_step,h_units=True))/(2.*z_step)*(1+z)
        else:
            # if z is between z_step/10 and z_step, reduce z_step to z, and then stick to two-sided derivative
            if (z > z_step/10.):
                z_step = z
                return (self.sigma(8,z-z_step,h_units=True)-self.sigma(8,z+z_step,h_units=True))/(2.*z_step)*(1+z)
            # if z is between 0 and z_step/10, use single-sided derivative with z_step/10
            else:
                z_step /=10
                return (self.sigma(8,z,h_units=True)-self.sigma(8,z+z_step,h_units=True))/z_step*(1+z)

    #################################
    # gives an estimation of f(z)*sigma8(z) at the scale of 8 h/Mpc, computed as (d sigma8/d ln a)
    def effective_f_sigma8_spline(self, z, Nz=20):
        """
        effective_f_sigma8_spline(z)

        Returns the time derivative of sigma8(z) computed as (d sigma8/d ln a)

        Parameters
        ----------
        z : float
                Desired redshift
        Nz : integer
                Number of values used to spline sigma8(z) in the range [z-0.1,z+0.1]

        Returns
        -------
        (d ln sigma8/d ln a)(z) (dimensionless)
        """

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)
        z_max = self.z_of_tau(np.exp(self.fo.ln_tau[0]))

        if (z<0) or (z>z_max):
            raise CosmoSevereError("You asked for effective_f_sigma8 at a redshift %e outside of the computed range [0,%e]"%(z,z_max))

        if (z<0.1):
            z_array = np.linspace(0, 0.2, num = Nz)
        elif (z<z_max-0.1):
            z_array = np.linspace(z-0.1, z+0.1, num = Nz)
        else:
            z_array = np.linspace(z_max-0.2, z_max, num = Nz)

        sig8_array = np.empty_like(z_array)
        for iz, zval in enumerate(z_array):
            sig8_array[iz] = self.sigma(8,zval,h_units=True)
        return -CubicSpline(z_array,sig8_array).derivative()(z)*(1+z)

   #################################
    def z_of_tau(self, tau):
        """
        Redshift corresponding to a given conformal time.

        Parameters
        ----------
        tau : float
                Conformal time
        """
        cdef double z
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        z = 1./pvecback[self.ba.index_bg_a]-1.

        free(pvecback)

        return z

    def Hubble(self, z):
        """
        Hubble(z)

        Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        H = pvecback[self.ba.index_bg_H]

        free(pvecback)

        return H

    def Om_m(self, z):
        """
        Omega_m(z)

        Return the matter density fraction (exactly, the quantity defined by Class as index_bg_Omega_m
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        Om_m = pvecback[self.ba.index_bg_Omega_m]

        free(pvecback)

        return Om_m

    def Om_b(self, z):
        """
        Omega_b(z)

        Return the baryon density fraction (exactly, the ratio of quantities defined by Class as
        index_bg_rho_b and index_bg_rho_crit in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        Om_b = pvecback[self.ba.index_bg_rho_b]/pvecback[self.ba.index_bg_rho_crit]

        free(pvecback)

        return Om_b

    def Om_cdm(self, z):
        """
        Omega_cdm(z)

        Return the cdm density fraction (exactly, the ratio of quantities defined by Class as
        index_bg_rho_cdm and index_bg_rho_crit in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        if self.ba.has_cdm == True:

            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            Om_cdm = pvecback[self.ba.index_bg_rho_cdm]/pvecback[self.ba.index_bg_rho_crit]

            free(pvecback)

        else:

            Om_cdm = 0.

        return Om_cdm

    def Om_ncdm(self, z):
        """
        Omega_ncdm(z)

        Return the ncdm density fraction (exactly, the ratio of quantities defined by Class as
        Sum_m [ index_bg_rho_ncdm1 + n ], with n=0...N_ncdm-1, and index_bg_rho_crit in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback

        if self.ba.has_ncdm == True:

            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            rho_ncdm = 0.
            for n in xrange(self.ba.N_ncdm):
                rho_ncdm += pvecback[self.ba.index_bg_rho_ncdm1+n]
            Om_ncdm = rho_ncdm/pvecback[self.ba.index_bg_rho_crit]

            free(pvecback)

        else:

            Om_ncdm = 0.

        return Om_ncdm

    def ionization_fraction(self, z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        xe = pvecthermo[self.th.index_th_xe]

        free(pvecback)
        free(pvecthermo)

        return xe

    def baryon_temperature(self, z):
        """
        baryon_temperature(z)

        Give the baryon temperature for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        Tb = pvecthermo[self.th.index_th_Tb]

        free(pvecback)
        free(pvecthermo)

        return Tb

    def T_cmb(self):
        """
        Return the CMB temperature
        """
        return self.ba.T_cmb

    # redundent with a previous Omega_m() funciton,
    # but we leave it not to break compatibility
    def Omega0_m(self):
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        return self.ba.Omega0_m

    def get_background(self):
        """
        Return an array of the background quantities at all times.

        Parameters
        ----------

        Returns
        -------
        background : dictionary containing background.
        """
        cdef char *titles
        cdef double* data
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if background_output_titles(&self.ba, titles)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.ba.bt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if background_output_data(&self.ba, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        background = {}

        for i in range(number_of_titles):
            background[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                background[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return background

    def get_thermodynamics(self):
        """
        Return the thermodynamics quantities.

        Returns
        -------
        thermodynamics : dictionary containing thermodynamics.
        """
        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if thermodynamics_output_titles(&self.ba, &self.th, titles)==_FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.th.tt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if thermodynamics_output_data(&self.ba, &self.th, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        thermodynamics = {}

        for i in range(number_of_titles):
            thermodynamics[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                thermodynamics[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return thermodynamics

    def get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """
        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if primordial_output_titles(&self.pt, &self.pm, titles)==_FAILURE_:
            raise CosmoSevereError(self.pm.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pm.lnk_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if primordial_output_data(&self.pt, &self.pm, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pm.error_message)

        primordial = {}

        for i in range(number_of_titles):
            primordial[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                primordial[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return primordial

    def get_perturbations(self):
        """
        Return scalar, vector and/or tensor perturbations as arrays for requested
        k-values.

        .. note::

            you need to specify both 'k_output_values', and have some
            perturbations computed, for instance by setting 'output' to 'tCl'.

        Returns
        -------
        perturbations : dict of array of dicts
                perturbations['scalar'] is an array of length 'k_output_values' of
                dictionary containing scalar perturbations.
                Similar for perturbations['vector'] and perturbations['tensor'].
        """

        perturbations = {}

        if self.pt.k_output_values_num<1:
            return perturbations

        cdef:
            Py_ssize_t j
            Py_ssize_t i
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            list names
            list tmparray
            dict tmpdict
            double[:,::1] data_mv
            double ** thedata
            int * thesizes

        # Doing the exact same thing 3 times, for scalar, vector and tensor. Sorry
        # for copy-and-paste here, but I don't know what else to do.
        for mode in ['scalar','vector','tensor']:
            if mode=='scalar' and self.pt.has_scalars:
                thetitles = <bytes> self.pt.scalar_titles
                thedata = self.pt.scalar_perturbations_data
                thesizes = self.pt.size_scalar_perturbation_data
            elif mode=='vector' and self.pt.has_vectors:
                thetitles = <bytes> self.pt.vector_titles
                thedata = self.pt.vector_perturbations_data
                thesizes = self.pt.size_vector_perturbation_data
            elif mode=='tensor' and self.pt.has_tensors:
                thetitles = <bytes> self.pt.tensor_titles
                thedata = self.pt.tensor_perturbations_data
                thesizes = self.pt.size_tensor_perturbation_data
            else:
                continue
            thetitles = str(thetitles.decode())
            names = thetitles.split("\t")[:-1]
            number_of_titles = len(names)
            tmparray = []
            if number_of_titles != 0:
                for j in range(self.pt.k_output_values_num):
                    timesteps = thesizes[j]//number_of_titles
                    tmpdict={}
                    data_mv = <double[:timesteps,:number_of_titles]> thedata[j]
                    for i in range(number_of_titles):
                        tmpdict[names[i]] = np.asarray(data_mv[:,i])
                    tmparray.append(tmpdict)
            perturbations[mode] = tmparray

        return perturbations

    def get_transfer(self, z=0., output_format='class'):
        """
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'mTk' and/or 'vCTk' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that 0<z<z_pk.

        Parameters
        ----------
        z  : redshift (default = 0)
        output_format  : ('class' or 'camb') Format transfer functions according to
                         CLASS convention (default) or CAMB convention.

        Returns
        -------
        tk : dictionary containing transfer functions.
        """
        cdef char *titles
        cdef double* data
        cdef char ic_info[1024]
        cdef FileName ic_suffix
        cdef file_format outf

        if (not self.pt.has_density_transfers) and (not self.pt.has_velocity_transfers):
            return {}

        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        index_md = self.pt.index_md_scalars;
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if perturbations_output_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pt.k_size[index_md]

        size_ic_data = timesteps*number_of_titles;
        ic_num = self.pt.ic_size[index_md];

        data = <double*>malloc(sizeof(double)*size_ic_data*ic_num)

        if perturbations_output_data_at_z(&self.ba, &self.pt, outf, <double> z, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        transfers = {}

        for index_ic in range(ic_num):
            if perturbations_output_firstline_and_ic_suffix(&self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
                raise CosmoSevereError(self.pt.error_message)
            ic_key = <bytes> ic_suffix

            tmpdict = {}
            for i in range(number_of_titles):
                tmpdict[names[i]] = np.zeros(timesteps, dtype=np.double)
                for index in range(timesteps):
                    tmpdict[names[i]][index] = data[index_ic*size_ic_data+index*number_of_titles+i]

            if ic_num==1:
                transfers = tmpdict
            else:
                transfers[ic_key] = tmpdict

        free(titles)
        free(data)

        return transfers

    def get_current_derived_parameters(self, names):
        """
        get_current_derived_parameters(names)

        Return a dictionary containing an entry for all the names defined in the
        input list.

        Parameters
        ----------
        names : list
                Derived parameters that can be asked from Monte Python, or
                elsewhere.

        Returns
        -------
        derived : dict

        .. warning::

            This method used to take as an argument directly the data class from
            Monte Python. To maintain compatibility with this old feature, a
            check is performed to verify that names is indeed a list. If not, it
            returns a TypeError. The old version of this function, when asked
            with the new argument, will raise an AttributeError.

        """
        if type(names) != type([]):
            raise TypeError("Deprecated")

        derived = {}
        for name in names:
            if name == 'h':
                value = self.ba.h
            elif name == 'H0':
                value = self.ba.h*100
            elif name == 'Omega0_lambda' or name == 'Omega_Lambda':
                value = self.ba.Omega0_lambda
            elif name == 'Omega0_fld':
                value = self.ba.Omega0_fld
            elif name == 'age':
                value = self.ba.age
            elif name == 'conformal_age':
                value = self.ba.conformal_age
            elif name == 'm_ncdm_in_eV':
                value = self.ba.m_ncdm_in_eV[0]
            elif name == 'm_ncdm_tot':
                value = self.ba.Omega0_ncdm_tot*self.ba.h*self.ba.h*93.14
            elif name == 'Neff':
                value = self.ba.Neff
            elif name == 'Omega_m':
                value = self.ba.Omega0_m
            elif name == 'omega_m':
                value = self.ba.Omega0_m*self.ba.h**2
            elif name == 'xi_idr':
                value = self.ba.T_idr/self.ba.T_cmb
            elif name == 'N_dg':
                value = self.ba.Omega0_idr/self.ba.Omega0_g*8./7.*pow(11./4.,4./3.)
            elif name == 'Gamma_0_nadm':
                value = self.th.a_idm_dr*(4./3.)*(self.ba.h*self.ba.h*self.ba.Omega0_idr)
            elif name == 'a_dark':
                value = self.th.a_idm_dr
            elif name == 'tau_reio':
                value = self.th.tau_reio
            elif name == 'z_reio':
                value = self.th.z_reio
            elif name == 'z_rec':
                value = self.th.z_rec
            elif name == 'tau_rec':
                value = self.th.tau_rec
            elif name == 'rs_rec':
                value = self.th.rs_rec
            elif name == 'rs_rec_h':
                value = self.th.rs_rec*self.ba.h
            elif name == 'ds_rec':
                value = self.th.ds_rec
            elif name == 'ds_rec_h':
                value = self.th.ds_rec*self.ba.h
            elif name == 'ra_rec':
                value = self.th.da_rec*(1.+self.th.z_rec)
            elif name == 'ra_rec_h':
                value = self.th.da_rec*(1.+self.th.z_rec)*self.ba.h
            elif name == 'da_rec':
                value = self.th.da_rec
            elif name == 'da_rec_h':
                value = self.th.da_rec*self.ba.h
            elif name == 'rd_rec':
                value = self.th.rd_rec
            elif name == 'z_star':
                value = self.th.z_star
            elif name == 'tau_star':
                value = self.th.tau_star
            elif name == 'rs_star':
                value = self.th.rs_star
            elif name == 'ds_star':
                value = self.th.ds_star
            elif name == 'ra_star':
                value = self.th.ra_star
            elif name == 'da_star':
                value = self.th.da_star
            elif name == 'rd_star':
                value = self.th.rd_star
            elif name == 'z_d':
                value = self.th.z_d
            elif name == 'tau_d':
                value = self.th.tau_d
            elif name == 'ds_d':
                value = self.th.ds_d
            elif name == 'ds_d_h':
                value = self.th.ds_d*self.ba.h
            elif name == 'rs_d':
                value = self.th.rs_d
            elif name == 'rs_d_h':
                value = self.th.rs_d*self.ba.h
            elif name == '100*theta_s':
                value = 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)
            elif name == '100*theta_star':
                value = 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)
            elif name == 'YHe':
                value = self.th.YHe
            elif name == 'n_e':
                value = self.th.n_e
            elif name == 'A_s':
                value = self.pm.A_s
            elif name == 'ln10^{10}A_s':
                value = log(1.e10*self.pm.A_s)
            elif name == 'n_s':
                value = self.pm.n_s
            elif name == 'alpha_s':
                value = self.pm.alpha_s
            elif name == 'beta_s':
                value = self.pm.beta_s
            elif name == 'r':
                # This is at the pivot scale
                value = self.pm.r
            elif name == 'r_0002':
                # at k_pivot = 0.002/Mpc
                value = self.pm.r*(0.002/self.pm.k_pivot)**(
                    self.pm.n_t-self.pm.n_s-1+0.5*self.pm.alpha_s*log(
                        0.002/self.pm.k_pivot))
            elif name == 'n_t':
                value = self.pm.n_t
            elif name == 'alpha_t':
                value = self.pm.alpha_t
            elif name == 'V_0':
                value = self.pm.V0
            elif name == 'V_1':
                value = self.pm.V1
            elif name == 'V_2':
                value = self.pm.V2
            elif name == 'V_3':
                value = self.pm.V3
            elif name == 'V_4':
                value = self.pm.V4
            elif name == 'epsilon_V':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                value = eps1*((1.-eps1/3.+eps2/6.)/(1.-eps1/3.))**2
            elif name == 'eta_V':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                eps23 = 1./8.*(self.pm.r**2/8.+(self.pm.n_s-1.)*self.pm.r-8.*self.pm.alpha_s)
                value = (2.*eps1-eps2/2.-2./3.*eps1**2+5./6.*eps1*eps2-eps2**2/12.-eps23/6.)/(1.-eps1/3.)
            elif name == 'ksi_V^2':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                eps23 = 1./8.*(self.pm.r**2/8.+(self.pm.n_s-1.)*self.pm.r-8.*self.pm.alpha_s)
                value = 2.*(1.-eps1/3.+eps2/6.)*(2.*eps1**2-3./2.*eps1*eps2+eps23/4.)/(1.-eps1/3.)**2
            elif name == 'exp_m_2_tau_As':
                value = exp(-2.*self.th.tau_reio)*self.pm.A_s
            elif name == 'phi_min':
                value = self.pm.phi_min
            elif name == 'phi_max':
                value = self.pm.phi_max
            elif name == 'sigma8':
                value = self.fo.sigma8[self.fo.index_pk_m]
            elif name == 'sigma8_cb':
                value = self.fo.sigma8[self.fo.index_pk_cb]
            elif name == 'k_eq':
                value = self.ba.a_eq*self.ba.H_eq
            elif name == 'a_eq':
                value = self.ba.a_eq
            elif name == 'z_eq':
                value = 1./self.ba.a_eq-1.
            elif name == 'H_eq':
                value = self.ba.H_eq
            elif name == 'tau_eq':
                value = self.ba.tau_eq
            elif name == 'g_sd':
                value = self.sd.sd_parameter_table[0]
            elif name == 'y_sd':
                value = self.sd.sd_parameter_table[1]
            elif name == 'mu_sd':
                value = self.sd.sd_parameter_table[2]
            elif name == 'network_delta_chi_squared':
                value = self.pt.network_deltachisquared
            else:
                raise CosmoSevereError("%s was not recognized as a derived parameter" % name)
            derived[name] = value
        return derived

    def nonlinear_scale(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_scale(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl_cb = np.zeros(z_size,'float64')
        #cdef double *k_nl
        #k_nl = <double*> calloc(z_size,sizeof(double))
        for index_z in range(z_size):
            if fourier_k_nl_at_z(&self.ba,&self.fo,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return k_nl

    def nonlinear_scale_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """

make        nonlinear_scale_cb(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size

        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl_cb = np.zeros(z_size,'float64')
        #cdef double *k_nl
        #k_nl = <double*> calloc(z_size,sizeof(double))
        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use pk, rather than pk_cb."
                )
        for index_z in range(z_size):
            if fourier_k_nl_at_z(&self.ba,&self.fo,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.fo.error_message)

        return k_nl_cb

    def fourier_hmcode_sigma8(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigma8_at_z(&self.ba,&self.fo,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_8

    def fourier_hmcode_sigma8_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigma8_at_z(&self.ba,&self.fo,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_8_cb

    def fourier_hmcode_sigmadisp(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp

    def fourier_hmcode_sigmadisp_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp_cb

    def fourier_hmcode_sigmadisp100(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp100(z, z_size)

        Return sigma_disp_100 for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp100_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp_100

    def fourier_hmcode_sigmadisp100_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmadisp100(z, z_size)

        Return sigma_disp_100 for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmadisp100_at_z(&self.ba,&self.fo,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_disp_100_cb

    def fourier_hmcode_sigmaprime(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmaprime(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmaprime_at_z(&self.ba,&self.fo,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_prime

    def fourier_hmcode_sigmaprime_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        fourier_hmcode_sigmaprime(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if fourier_hmcode_sigmaprime_at_z(&self.ba,&self.fo,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.fo.error_message)

        return sigma_prime_cb

    def __call__(self, ctx):
        """
        Function to interface with CosmoHammer

        Parameters
        ----------
        ctx : context
                Contains several dictionaries storing data and cosmological
                information

        """
        data = ctx.get('data')  # recover data from the context

        # If the module has already been called once, clean-up
        if self.state:
            self.struct_cleanup()

        # Set the module to the current values
        self.set(data.cosmo_arguments)
        self.compute(["lensing"])

        # Compute the derived paramter value and store them
        params = ctx.getData()
        self.get_current_derived_parameters(
            data.get_mcmc_parameters(['derived']))
        for elem in data.get_mcmc_parameters(['derived']):
            data.mcmc_parameters[elem]['current'] /= \
                data.mcmc_parameters[elem]['scale']
            params[elem] = data.mcmc_parameters[elem]['current']

        ctx.add('boundary', True)
        # Store itself into the context, to be accessed by the likelihoods
        ctx.add('cosmo', self)

    def get_pk_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk

    def get_pk_cb_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk_cb

    def Omega0_k(self):
        """ Curvature contribution """
        return self.ba.Omega0_k

    def Omega0_cdm(self):
        return self.ba.Omega0_cdm

    def spectral_distortion_amplitudes(self):
        if self.sd.type_size == 0:
          raise CosmoSevereError("No spectral distortions have been calculated. Check that the output contains 'Sd' and the compute level is at least 'distortions'.")
        cdef np.ndarray[DTYPE_t, ndim=1] sd_type_amps = np.zeros(self.sd.type_size,'float64')
        for i in range(self.sd.type_size):
          sd_type_amps[i] = self.sd.sd_parameter_table[i]
        return sd_type_amps

    def spectral_distortion(self):
        if self.sd.x_size == 0:
          raise CosmoSevereError("No spectral distortions have been calculated. Check that the output contains 'Sd' and the compute level is at least 'distortions'.")
        cdef np.ndarray[DTYPE_t, ndim=1] sd_amp = np.zeros(self.sd.x_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sd_nu = np.zeros(self.sd.x_size,'float64')
        for i in range(self.sd.x_size):
          sd_amp[i] = self.sd.DI[i]*self.sd.DI_units*1.e26
          sd_nu[i] = self.sd.x[i]*self.sd.x_to_nu
        return sd_nu,sd_amp

    ################################################################################
    # utility functions for neural networks
    ################################################################################
    def get_tau_source(self):
        """
        Return the tau array on which the source functions are sampled.

        Returns
        -------
        tau: numpy array
        """
        cdef:
            int i_tau;
            double * tau = self.pt.tau_sampling;
            int tau_size = self.pt.tau_size
            double [:] numpy_tau = np.zeros(tau_size, dtype=np.double)

        for i_tau in range(tau_size):
            numpy_tau[i_tau] = tau[i_tau]

        return np.asarray(numpy_tau)

    def get_k_tau(self):
        """
        Return the k, tau grid values of the source functions.

        Returns
        -------
        k_array : numpy array containing k values.
        tau_array: numpy array containing tau values.
        """
        cdef:
            int i_k;
            int i_tau;
            int index_md = self.pt.index_md_scalars;
            double * k = self.pt.k[index_md];
            double * tau = self.pt.tau_sampling;
            int k_size = self.pt.k_size[index_md];
            int tau_size = self.pt.tau_size
            double [:] numpy_k = np.zeros(k_size,dtype=np.double)
            double [:] numpy_tau = np.zeros(tau_size,dtype=np.double)

        for i_k in range(k_size):
            numpy_k[i_k] = k[i_k]
        for i_tau in range(tau_size):
            numpy_tau[i_tau] = tau[i_tau]

        return np.asarray(numpy_k), np.asarray(numpy_tau)

    def get_quantities_at_RM_equality(self):
        return self.ba.tau_eq, self.ba.a_eq, self.ba.H_eq

    def get_bg_z(self):
        cdef:
            int index_tau;
            int bt_size = self.ba.bt_size;
            double * z_table = self.ba.z_table;
            double [:] np_z_table = np.zeros((bt_size));

        for index_tau in range(bt_size):
            np_z_table[index_tau] = z_table[index_tau]

        return np.asarray(np_z_table)

    def get_bg_tau(self):
        cdef:
            int index_tau;
            int bt_size = self.ba.bt_size;
            double * tau_table = self.ba.tau_table;
            double [:] np_tau_table = np.zeros((bt_size));

        for index_tau in range(bt_size):
            np_tau_table[index_tau] = tau_table[index_tau]

        return np.asarray(np_tau_table)

    def get_backgrounds_for_NN(self):
        """
        Return the comoving sound horizon of photons from the background module and the
        corresponding conformal time values. We also need to calculate Omega_m
        """
        cdef:
            int index_tau;
            int bt_size = self.ba.bt_size;
            int bg_size = self.ba.bg_size;
            int index_bg_rs = self.ba.index_bg_rs;
            double * background_table = self.ba.background_table;
            double [:] r_s = np.empty(bt_size, dtype=np.double);
            double [:] rho_b = np.empty(bt_size, dtype=np.double);
            double [:] rho_g = np.empty(bt_size, dtype=np.double);
            double [:] tau_bg = np.empty(bt_size, dtype=np.double);
            double [:] a = np.empty(bt_size, dtype=np.double);
            double [:] H = np.empty(bt_size, dtype=np.double);
            double [:] D = np.empty(bt_size, dtype=np.double);
            double [:] Omega_m = np.empty(bt_size, dtype=np.double);


        for index_tau in range(bt_size):
            tau_bg[index_tau] = self.ba.tau_table[index_tau]

            r_s[index_tau] = background_table[index_tau*bg_size+index_bg_rs]
            rho_b[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_rho_b];
            rho_g[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_rho_g];
            a[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_a];
            H[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_H];
            D[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_D];
            Omega_m[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_Omega_m];
            



        return {
                "tau": np.asarray(tau_bg),
                "r_s": np.asarray(r_s),
                "rho_b": np.asarray(rho_b),
                "rho_g": np.asarray(rho_g),
                "a": np.asarray(a),
                "H": np.asarray(H),
                "D": np.asarray(D),
                "Omega_m": np.asarray(Omega_m),
                }

    def tau_of_z(self,z):
        cdef double tau
        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)
        return tau

    def get_z_split_eisw_lisw(self):
        return self.pt.eisw_lisw_split_z

    def get_tau_split_eisw_lisw(self):
        return self.tau_of_z(self.get_z_split_eisw_lisw())

    @cython.boundscheck(False)
    def get_thermos_for_NN(self):
        """
        Return the photon comoving damping scale, visibility function, its conformal time
        derivative and the corresponding conformal time values.
        """
        cdef:
            double tau;
            double * z_table = self.th.z_table;
            double *tau_table = self.th.tau_table;
            double * thermodynamics_table = self.th.thermodynamics_table;
            int th_size = self.th.th_size;
            int index_z;
            int tt_size = self.th.tt_size;
            int index_th_r_d = self.th.index_th_r_d;
            int index_th_g = self.th.index_th_g;
            int index_th_dg = self.th.index_th_dg;
            int index_th_dg_reco = self.th.index_th_dg_reco;
            int index_th_dg_reio = self.th.index_th_dg_reio;
            int index_th_dkappa = self.th.index_th_dkappa;
            int index_th_exp_m_kappa = self.th.index_th_exp_m_kappa;
            double [:] numpy_r_d = np.empty(tt_size, dtype=np.double);
            double [:] numpy_g = np.empty(tt_size, dtype=np.double);
            double [:] numpy_g_reco = np.empty(tt_size, dtype=np.double);
            double [:] numpy_g_reio = np.empty(tt_size, dtype=np.double);
            double [:] numpy_dg = np.empty(tt_size, dtype=np.double);
            double [:] numpy_dg_reco = np.empty(tt_size, dtype=np.double);
            double [:] numpy_dg_reio = np.empty(tt_size, dtype=np.double);
            double [:] numpy_e_kappa = np.empty(tt_size, dtype=np.double);
            double [:] numpy_dkappa = np.empty(tt_size, dtype=np.double);
            double [:] numpy_tau = np.empty(tt_size, dtype=np.double);     #this takes 2.5e-5 sec

        for index_z in prange(tt_size, nogil=True):

            # if background_tau_of_z(&self.ba,z_table[index_z],&tau)==_FAILURE_:
            #     raise CosmoSevereError(self.ba.error_message)

            numpy_tau[index_z] = tau_table[index_z]
            numpy_r_d[index_z] = thermodynamics_table[index_z*th_size + index_th_r_d]
            numpy_g[index_z] = thermodynamics_table[index_z*th_size + index_th_g]
            numpy_g_reco[index_z] = thermodynamics_table[index_z*th_size + self.th.index_th_g_reco]
            numpy_g_reio[index_z] = thermodynamics_table[index_z*th_size + self.th.index_th_g_reio]
            numpy_dg[index_z] = thermodynamics_table[index_z*th_size + index_th_dg]
            numpy_dg_reco[index_z] = thermodynamics_table[index_z*th_size + index_th_dg_reco]
            numpy_dg_reio[index_z] = thermodynamics_table[index_z*th_size + index_th_dg_reio]
            numpy_dkappa[index_z] = thermodynamics_table[index_z*th_size + index_th_dkappa]
            numpy_e_kappa[index_z] = thermodynamics_table[index_z*th_size + index_th_exp_m_kappa]     #this takes 4e-4

        g_reco = np.asarray(numpy_g_reco)
        g_reio = np.asarray(numpy_g_reio)
        g_reco_prime = np.asarray(numpy_dg_reco)
        g_reio_prime = np.asarray(numpy_dg_reio)

        ret = {
                "r_d": np.asarray(numpy_r_d),
                "g": np.asarray(numpy_g),
                "g_prime": np.asarray(numpy_dg),
                "g_reco": g_reco,
                "g_reco_prime": g_reco_prime,
                "g_reio": g_reio,
                "g_reio_prime": g_reio_prime,
                "e_kappa": np.asarray(numpy_e_kappa),
                "tau": np.asarray(numpy_tau),
                "dkappa": np.asarray(numpy_dkappa),
                } #takes 1e-5 sec

        return ret

    def a_of_tau(self, tau):
        """
        a(tau)

        Return the scale factor

        Parameters
        ----------
        tau : float
                Desired conformal time
        """
        cdef int last_index
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size_short,sizeof(double))
        if background_at_tau(&self.ba,tau,short_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        a = pvecback[self.ba.index_bg_a]

        free(pvecback)

        return a


    def get_sources(self):
        """
        Return the source functions for all k, tau in the grid.

        Returns
        -------
        sources : dictionary containing source functions.
        k_array : numpy array containing k values.
        tau_array: numpy array containing tau values.
        """
        sources = {}

        cdef:
            int index_k, index_tau, i_index_type;
            int index_type;
            int index_md = self.pt.index_md_scalars;
            double * k = self.pt.k[index_md];
            double * tau = self.pt.tau_sampling;
            int index_ic = self.pt.index_ic_ad;
            int k_size = self.pt.k_size[index_md];
            int tau_size = self.pt.tau_size;
            int tp_size = self.pt.tp_size[index_md]
            double *** sources_ptr = self.pt.sources
            double data
            double [:,:] tmparray = np.zeros((k_size, tau_size))
            double [:] k_array = np.zeros(k_size)
            double [:] tau_array = np.zeros(tau_size)

        names = []

        for index_k in range(k_size):
            k_array[index_k] = k[index_k]
        for index_tau in range(tau_size):
            tau_array[index_tau] = tau[index_tau]

        # indices = []

        # if self.pt.has_source_t:
        #     indices.extend([
        #         self.pt.index_tp_t0,
        #         self.pt.index_tp_t1,
        #         self.pt.index_tp_t2
        #         ])
        #     names.extend([
        #         "t0",
        #         "t1",
        #         "t2"
        #         ])
        # if self.pt.has_source_p:
        #     indices.append(self.pt.index_tp_p)
        #     names.append("p")
        # if self.pt.has_source_phi:
        #     indices.append(self.pt.index_tp_phi)
        #     names.append("phi")
        # if self.pt.has_source_phi_plus_psi:
        #     indices.append(self.pt.index_tp_phi_plus_psi)
        #     names.append("phi_plus_psi")
        # if self.pt.has_source_phi_prime:
        #     indices.append(self.pt.index_tp_phi_prime)
        #     names.append("phi_prime")
        # if self.pt.has_source_psi:
        #     indices.append(self.pt.index_tp_psi)
        #     names.append("psi")
        # if self.pt.has_source_H_T_Nb_prime:
        #     indices.append(self.pt.index_tp_H_T_Nb_prime)
        #     names.append("H_T_Nb_prime")
        # if self.pt.index_tp_k2gamma_Nb:
        #     indices.append(self.pt.index_tp_k2gamma_Nb)
        #     names.append("k2gamma_Nb")
        # if self.pt.has_source_h:
        #     indices.append(self.pt.index_tp_h)
        #     names.append("h")
        # if self.pt.has_source_h_prime:
        #     indices.append(self.pt.index_tp_h_prime)
        #     names.append("h_prime")
        # if self.pt.has_source_eta:
        #     indices.append(self.pt.index_tp_eta)
        #     names.append("eta")
        # if self.pt.has_source_eta_prime:
        #     indices.append(self.pt.index_tp_eta_prime)
        #     names.append("eta_prime")
        # if self.pt.has_source_delta_tot:
        #     indices.append(self.pt.index_tp_delta_tot)
        #     names.append("delta_tot")
        # if self.pt.has_source_delta_m:
        #     indices.append(self.pt.index_tp_delta_m)
        #     names.append("delta_m")
        # if self.pt.has_source_delta_cb:
        #     indices.append(self.pt.index_tp_delta_cb)
        #     names.append("delta_cb")
        # if self.pt.has_source_delta_g:
        #     indices.append(self.pt.index_tp_delta_g)
        #     names.append("delta_g")
        # if self.pt.has_source_delta_b:
        #     indices.append(self.pt.index_tp_delta_b)
        #     names.append("delta_b")
        # if self.pt.has_source_delta_cdm:
        #     indices.append(self.pt.index_tp_delta_cdm)
        #     names.append("delta_cdm")
        # if self.pt.has_source_delta_idm:
        #     indices.append(self.pt.index_tp_delta_idm)
        #     names.append("delta_idm")
        # if self.pt.has_source_delta_dcdm:
        #     indices.append(self.pt.index_tp_delta_dcdm)
        #     names.append("delta_dcdm")
        # if self.pt.has_source_delta_fld:
        #     indices.append(self.pt.index_tp_delta_fld)
        #     names.append("delta_fld")
        # if self.pt.has_source_delta_scf:
        #     indices.append(self.pt.index_tp_delta_scf)
        #     names.append("delta_scf")
        # if self.pt.has_source_delta_dr:
        #     indices.append(self.pt.index_tp_delta_dr)
        #     names.append("delta_dr")
        # if self.pt.has_source_delta_ur:
        #     indices.append(self.pt.index_tp_delta_ur)
        #     names.append("delta_ur")
        # if self.pt.has_source_delta_idr:
        #     indices.append(self.pt.index_tp_delta_idr)
        #     names.append("delta_idr")
        # if self.pt.has_source_delta_ncdm:
        #     for incdm in range(self.ba.N_ncdm):
        #       indices.append(self.pt.index_tp_delta_ncdm1+incdm)
        #       names.append("delta_ncdm[{}]".format(incdm))
        # if self.pt.has_source_theta_tot:
        #     indices.append(self.pt.index_tp_theta_tot)
        #     names.append("theta_tot")
        # if self.pt.has_source_theta_m:
        #     indices.append(self.pt.index_tp_theta_m)
        #     names.append("theta_m")
        # if self.pt.has_source_theta_cb:
        #     indices.append(self.pt.index_tp_theta_cb)
        #     names.append("theta_cb")
        # if self.pt.has_source_theta_g:
        #     indices.append(self.pt.index_tp_theta_g)
        #     names.append("theta_g")
        # if self.pt.has_source_theta_b:
        #     indices.append(self.pt.index_tp_theta_b)
        #     names.append("theta_b")
        # if self.pt.has_source_theta_cdm:
        #     indices.append(self.pt.index_tp_theta_cdm)
        #     names.append("theta_cdm")
        # if self.pt.has_source_theta_idm:
        #     indices.append(self.pt.index_tp_theta_idm)
        #     names.append("theta_idm")
        # if self.pt.has_source_theta_dcdm:
        #     indices.append(self.pt.index_tp_theta_dcdm)
        #     names.append("theta_dcdm")
        # if self.pt.has_source_theta_fld:
        #     indices.append(self.pt.index_tp_theta_fld)
        #     names.append("theta_fld")
        # if self.pt.has_source_theta_scf:
        #     indices.append(self.pt.index_tp_theta_scf)
        #     names.append("theta_scf")
        # if self.pt.has_source_theta_dr:
        #     indices.append(self.pt.index_tp_theta_dr)
        #     names.append("theta_dr")
        # if self.pt.has_source_theta_ur:
        #     indices.append(self.pt.index_tp_theta_ur)
        #     names.append("theta_ur")
        # if self.pt.has_source_theta_idr:
        #     indices.append(self.pt.index_tp_theta_idr)
        #     names.append("theta_idr")
        # if self.pt.has_source_theta_ncdm:
        #     for incdm in range(self.ba.N_ncdm):
        #       indices.append(self.pt.index_tp_theta_ncdm1+incdm)
        #       names.append("theta_ncdm[{}]".format(incdm))

        # for index_type, name in zip(indices, names):
        #     tmparray = np.empty((k_size,tau_size))

        indices = {}

        if self.pt.has_source_t:
            indices.update({
                "t0": self.pt.index_tp_t0,
                # "t0_sw": self.pt.index_tp_t0_sw,
                "t0_isw": self.pt.index_tp_t0_isw,
                # "t0_reco": self.pt.index_tp_t0_reco,
                # "t0_reio": self.pt.index_tp_t0_reio,
                "t0_reco_no_isw": self.pt.index_tp_t0_reco_no_isw,
                "t0_reio_no_isw": self.pt.index_tp_t0_reio_no_isw,
                "t1": self.pt.index_tp_t1,
                "t2": self.pt.index_tp_t2,
                "t2_reco": self.pt.index_tp_t2_reco,
                "t2_reio": self.pt.index_tp_t2_reio
            })
        if self.pt.has_source_p:
            indices["p"] = self.pt.index_tp_p
        if self.pt.has_source_delta_m:
            indices["delta_m"] = self.pt.index_tp_delta_m
        if self.pt.has_source_delta_cb:
            indices["delta_cb"] = self.pt.index_tp_delta_cb
        if self.pt.has_source_delta_g:
            indices["delta_g"] = self.pt.index_tp_delta_g
        if self.pt.has_source_theta_m:
            indices["theta_m"] = self.pt.index_tp_theta_m
        if self.pt.has_source_theta_b:
            indices["theta_b"] = self.pt.index_tp_theta_b
        if self.pt.has_source_phi:
            indices["phi"] = self.pt.index_tp_phi
        if self.pt.has_source_phi_plus_psi:
            indices["phi_plus_psi"] = self.pt.index_tp_phi_plus_psi
        if self.pt.has_source_phi_prime:
            indices["phi_prime"] = self.pt.index_tp_phi_prime
        if self.pt.has_source_psi:
            indices["psi"] = self.pt.index_tp_psi

        for name, index_type in indices.items():
            for index_k in range(k_size):
                for index_tau in range(tau_size):
                    tmparray[index_k][index_tau] = sources_ptr[index_md][index_ic*tp_size+index_type][index_tau*k_size + index_k];

            sources[name] = np.asarray(tmparray)
            tmparray = np.zeros((k_size,tau_size))

        return (sources, np.asarray(k_array), np.asarray(tau_array))

    def nn_workspace(self):
        # TODO don't do this here
        workspace = self._pars["workspace_path"]
        if any(isinstance(workspace, t) for t in [str, bytes, os.PathLike]):
            # Check whether any generations are requested for nn
            joint = list(set(self._pars).intersection(classynet.models.ALL_NETWORK_STRINGS))
            if len(joint)>0:
                self.NN_generations = {name : self._pars[name] for name in joint}
                workspace = classynet.workspace.GenerationalWorkspace(workspace,self.NN_generations)
            else:
                workspace = classynet.workspace.Workspace(workspace)
        return workspace

    def nn_debug_enabled(self):
        return bool(self._pars.get("nn_debug", False))

    def k_min(self):
        """
        k_min as determined by K (taken from perturbations.c:perturb_get_k_list)
        """
        if self.ba.sgnK == 0:
            # K<0 (flat)  : start close to zero
            return self.pr.k_min_tau0 / self.ba.conformal_age
        elif self.ba.sgnK == -1:
            # K<0 (open)  : start close to sqrt(-K)
            # (in transfer modules, for scalars, this will correspond to q close to zero;
            # for vectors and tensors, this value is even smaller than the minimum necessary value)
            return np.sqrt(-self.ba.K + pow(self.pr.k_min_tau0 / self.ba.conformal_age / self.th.angular_rescaling, 2))
        elif self.ba.sgnK == 1:
            # K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K)
            return np.sqrt((8.-1.e-4) * self.ba.K);
        else:
            raise ValueError("Unrecognized value of K = {}!".format(self.ba.K))

    def get_q(self):
        """
        Get the q list from transfer module.
        """
        cdef:
            int i_q
            double * q = self.tr.q
            int q_size = self.tr.q_size
            double [:] numpy_q = np.empty(q_size, dtype=np.double)

        for i_q in range(q_size):
            numpy_q[i_q] = q[i_q]

        return np.asarray(numpy_q)

    def translate_source_to_index(self, name):
        mapping = {
            "t0":             self.pt.index_tp_t0,
            "t1":             self.pt.index_tp_t1,
            "t2":             self.pt.index_tp_t2,
            "t0_reco_no_isw": self.pt.index_tp_t0_reco_no_isw,
            "t0_reio_no_isw": self.pt.index_tp_t0_reio_no_isw,
            "t0_isw":         self.pt.index_tp_t0_isw,
            "t2_reco":        self.pt.index_tp_t2_reco,
            "t2_reio":        self.pt.index_tp_t2_reio,
            "phi_plus_psi":   self.pt.index_tp_phi_plus_psi,
            "delta_m":        self.pt.index_tp_delta_m,
            "delta_cb":       self.pt.index_tp_delta_cb,
            "p":              self.pt.index_tp_p,
        }
        return mapping[name]
