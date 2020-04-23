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
from math import exp,log
import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i

# Import the .pxd containing definitions
from cclassy cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

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
    cdef thermo th
    cdef perturbs pt
    cdef primordial pm
    cdef nonlinear nl
    cdef transfers tr
    cdef spectra sp
    cdef output op
    cdef lensing le
    cdef distortions sd
    cdef file_content fc

    cpdef int ready # Flag to see if classy can currently compute
    cpdef int allocated # Flag to see if classy structs are allocated already
    cpdef object _pars # Dictionary of the parameters
    cpdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

    cpdef bint use_NN
    cpdef object predictor

    # Defining two new properties to recover, respectively, the parameters used
    # or the age (set after computation). Follow this syntax if you want to
    # access other quantities. Alternatively, you can also define a method, and
    # call it (see _T_cmb method, at the very bottom).
    property pars:
        def __get__(self):
            return self._pars
    property state:
        def __get__(self):
            return self.ready
    property Omega_nu:
        def __get__(self):
            return self.ba.Omega0_ncdm_tot
    property nonlinear_method:
        def __get__(self):
            return self.nl.method

    def set_default(self):
        _pars = {
            "output":"tCl mPk",}
        self.set(**_pars)

    def __cinit__(self, default=False):
        cpdef char* dumc
        self.ready = False
        self.allocated = False
        self.use_NN = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()

    def enable_NN(self, predictor):
        """
        This function must be called if the user wants to use NNs to predict
        the source functions. In that case, a Predictor instance must be passed
        which will predict the source functions.
        """
        self.use_NN = True
        self.predictor = predictor
        # If using neural networks, set 'bad' parameters
        # (explanation in method definition) before calling
        # compute
        self._set_parameters_for_NN()

    # Set up the dictionary
    def set(self,*pars,**kars):
        if len(pars)==1:
            self._pars.update(dict(pars[0]))
        elif len(pars)!=0:
            raise CosmoSevereError("bad call")

        self._pars.update(kars)
        self.ready=False
        return True

    def empty(self):
        self._pars = {}
        self.ready = False
        self.use_NN = False
        self.predictor = None

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
            dumcp = str(self._pars[kk]).encode()
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = _FALSE_
            i+=1

    # Called at the end of a run, to free memory
    def struct_cleanup(self):
        if self.ready == _FALSE_:
             return
        if "distortions" in self.ncp:
            distortions_free(&self.sd)
        if "lensing" in self.ncp:
            lensing_free(&self.le)
        if "spectra" in self.ncp:
            spectra_free(&self.sp)
        if "transfer" in self.ncp:
            transfer_free(&self.tr)
        if "nonlinear" in self.ncp:
            nonlinear_free(&self.nl)
        if "primordial" in self.ncp:
            primordial_free(&self.pm)
        if "perturb" in self.ncp:
            perturb_free(&self.pt)
        if "thermodynamics" in self.ncp:
            thermodynamics_free(&self.th)
        if "background" in self.ncp:
            background_free(&self.ba)
        self.ready = False
        self.allocated = False

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
            level.append("lensing")
        if "lensing" in level:
            level.append("spectra")
        if "spectra" in level:
            level.append("transfer")
        if "transfer" in level:
            level.append("nonlinear")
        if "nonlinear" in level:
            level.append("primordial")
        if "primordial" in level:
            level.append("perturb")
        if "perturb" in level:
            level.append("thermodynamics")
        if "thermodynamics" in level:
            level.append("background")
        if len(level)!=0 :
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

    def _set_parameters_for_NN(self):
        """
        This method will perform two things:

        1.
        It will set 'bad' parameters for the perturbation module in order to
        'bypass' the full calculation but still initialize the structs, etc.
        This needs to be done (TODO: for now...) when using NNs to predict
        source functions.

        2.
        It will set additional parameters (see below as `additional_params`) that
        are required for NN evaluation.
        """

        bad_params = {
            "l_max_g": 4,
            "l_max_ur": 4,
            "l_max_pol_g": 4,
            "l_max_ncdm": 4,
            # "l_max_dr": 4,
            # "l_max_g_ten": 4,
            # "l_max_pol_g_ten": 4,

            "perturb_integration_stepsize": 5,
            "tol_perturb_integration": 1e10,
            # "perturb_sampling_stepsize": 100,
            "radiation_streaming_approximation": 0,
            "radiation_streaming_trigger_tau_over_tau_k": 0.01,
            "radiation_streaming_trigger_tau_c_over_tau": 0.2,
            "ur_fluid_approximation": 0,
            "ur_fluid_trigger_tau_over_tau_k": 0.3,
            # "neglect_CMB_sources_below_visibility": 1e5
        }
        self.set(bad_params)

    cdef void overwrite_source_function(self, int index_md, int index_ic,
            int index_type,
            int k_NN_size, int tau_size, double[:, :] S):
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
        # Allocate memory for source function
        self.pt.sources[index_md][index_ic * tp_size + index_type] = <double*> malloc(k_NN_size * tau_size * sizeof(double))

        for index_tau in range(tau_size):
            for index_k in range(k_NN_size):
                # GS TODO: More efficient way to copy memory from numpy array??
                self.pt.sources[index_md][index_ic*tp_size + index_type][index_tau*k_NN_size + index_k] = S[index_k][index_tau]

    def compute(self, level=["distortions"], performance_report=None):
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

        timer = Timer()

        timer.start("compute")

        # Append to the list level all the modules necessary to compute.
        level = self._check_task_dependency(level)

        # Check if this function ran before (self.ready should be true), and
        # if no other modules were requested, i.e. if self.ncp contains (or is
        # equivalent to) level. If it is the case, simply stop the execution of
        # the function.
        if self.ready and self.ncp.issuperset(level):
            return

        # Check if already allocated to prevent memory leaks
        if self.allocated:
            self.struct_cleanup()

        # Otherwise, proceed with the normal computation.
        self.ready = False

        # Equivalent of writing a parameter file
        self._fillparfile()

        # self.ncp will contain the list of computed modules (under the form of
        # a set, instead of a python list)
        self.ncp=set()

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
                                    &self.pt, &self.tr, &self.pm, &self.sp,
                                    &self.nl, &self.le, &self.sd, &self.op, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.ncp.add("input")
            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            for i in range(self.fc.size):
                if self.fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(self.fc.name[i].decode())
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
            int [:] index_types = np.zeros(tot_num_of_sources,dtype=np.int32)
            int [:] index_types_mask = np.zeros((tot_num_of_sources),dtype=np.int32)
            # double [:,:] NN_interpolated
            # double [:] NN_interpolated
            # TODO remove some of the unused ones here
            double [:, :, :] NN_prediction
            double * c_NN_sources

        if "perturb" in level:

            timer.start("perturb")
            timer.start("perturb_init")


            # Allocate memory for ALL source functions (since transfer.c iterates over them)

            if self.use_NN:
              self.pt.perform_NN_skip = _TRUE_
            if perturb_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            self.ncp.add("perturb")
            timer.end("perturb_init")

            # flag for using NN
            if self.use_NN:
                timer.start("neural network complete")
                timer.start("neural network initialization")

                index_md = self.pt.index_md_scalars;
                k_size = self.pt.k_size[index_md];
                tau_size = self.pt.tau_size;
                # tau_NN_size = len(self.tau_NN);
                index_ic = self.pt.index_ic_ad;

                tp_size = self.pt.tp_size[index_md];


                # c_tau_NN = <double*>malloc(tau_NN_size*sizeof(double))


                tau_CLASS = np.zeros((tau_size))
                for index_tau in range(tau_size):
                    tau_CLASS[index_tau] = self.pt.tau_sampling[index_tau]

                requested_index_types = []

                # add all sources that class calculates to the list for predicting
                if self.pt.has_source_t:
                    index_types[0] = self.pt.index_tp_t0
                    index_types[1] = self.pt.index_tp_t1
                    index_types[2] = self.pt.index_tp_t2
                    index_types_mask[0:3] = 1
                    requested_index_types.extend([
                        self.pt.index_tp_t0,
                        self.pt.index_tp_t1,
                        self.pt.index_tp_t2,
                        ])
                if self.pt.has_source_p:
                    index_types[3] = self.pt.index_tp_p
                    index_types_mask[3] = 1
                    # TODO Explicitly add self.pt.index_tp_p to requested_index_types
                    # or compute from ST2 here later?
                if self.pt.has_source_phi_plus_psi:
                    index_types[4] = self.pt.index_tp_phi_plus_psi
                    index_types_mask[4] = 1
                    requested_index_types.append(self.pt.index_tp_phi_plus_psi)

                if self.pt.has_source_delta_m:
                    index_types[5] = self.pt.index_tp_delta_m
                    index_types_mask[5] = 1

                    requested_index_types.append(self.pt.index_tp_delta_m)

                '''
                if self.pt.has_source_delta_g:
                    index_type_list.append(self.pt.index_tp_delta_g)
                    names.append('delta_g')
                if self.pt.has_source_theta_m:
                    index_type_list.append(self.pt.index_tp_theta_m)
                    names.append('theta_m')
                if self.pt.has_source_phi:
                    index_type_list.append(self.pt.index_tp_phi)
                    names.append('phi')
                if self.pt.has_source_phi_prime:
                    index_type_list.append(self.pt.index_tp_phi_prime)
                    names.append('phi_prime')
                if self.pt.has_source_psi:
                    index_type_list.append(self.pt.index_tp_psi)
                    names.append('psi')
                '''

                timer.end("neural network initialization")

                timer.start("get all sources")
                k_NN, NN_prediction = self.predictor.predict_all(self, np.asarray(tau_CLASS), add_k0=True)
                timer.end("get all sources")

                timer.start("overwrite k array")

                # Copy k values from NN
                k_NN_size = len(k_NN)
                free(self.pt.k[index_md])
                self.pt.k[index_md] = <double*>malloc(k_NN_size * sizeof(double))
                for i_k in range(k_NN_size):
                    self.pt.k[index_md][i_k] = k_NN[i_k]
                self.pt.k_size[index_md] = k_NN_size
                self.pt.k_size_cl[index_md] = k_NN_size

                timer.end("overwrite k array")

                timer.start("allocate unused source functions")
                for index_type in range(tp_size):
                    # Using malloc instead of calloc here will cause the splining
                    # in transfer.c to explode, but that doesn't seem to be an issue.
                    # Using malloc over calloc saves about a factor of 10 in runtime.
                    # self.pt.sources[index_md][index_ic * tp_size + index_type] = <double*> calloc(k_NN_size * tau_size,  sizeof(double))
                    self.pt.sources[index_md][index_ic * tp_size + index_type] = <double*> malloc(k_NN_size * tau_size * sizeof(double))
                timer.end("allocate unused source functions")

                timer["neural network evaluation"] = self.predictor.time_prediction
                timer["neural network input transformation"] = self.predictor.time_input_transformation
                timer["neural network output transformation"] = self.predictor.time_output_transformation

                for key, value in self.predictor.time_prediction_per_network.items():
                    timer["indiv. network: '{}'".format(key)] = value

                timer.start("overwrite source functions")

                c_NN_sources = <double*>malloc(k_NN_size * tau_size * sizeof(double))


                ############################################################
                if self.pt.has_source_t:
                    self.overwrite_source_function(
                            index_md, index_ic,
                            self.pt.index_tp_t0,
                            k_NN_size, tau_size, NN_prediction[0, :, :]
                            )
                    self.overwrite_source_function(
                            index_md, index_ic,
                            self.pt.index_tp_t1,
                            k_NN_size, tau_size, NN_prediction[1, :, :]
                            )
                    self.overwrite_source_function(
                            index_md, index_ic,
                            self.pt.index_tp_t2,
                            k_NN_size, tau_size, NN_prediction[2, :, :]
                            )

                if self.pt.has_source_p:
                    self.overwrite_source_function(
                            index_md, index_ic,
                            self.pt.index_tp_p,
                            k_NN_size, tau_size, np.sqrt(6) * NN_prediction[2, :, :]
                            )

                if self.pt.has_source_phi_plus_psi:
                    self.overwrite_source_function(
                            index_md, index_ic,
                            self.pt.index_tp_phi_plus_psi,
                            k_NN_size, tau_size, NN_prediction[3, :, :]
                            )

                if self.pt.has_source_delta_m:
                    self.overwrite_source_function(
                            index_md, index_ic,
                            self.pt.index_tp_delta_m,
                            k_NN_size, tau_size, NN_prediction[4, :, :]
                            )


                free(c_NN_sources)
                timer.end("overwrite source functions")
                ############################################################

                timer.end("neural network complete")

            timer.end("perturb")

        if "primordial" in level:
            timer.start("primordial")
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            self.ncp.add("primordial")
            timer.end("primordial")

        if "nonlinear" in level:
            timer.start("nonlinear")
            if nonlinear_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.nl) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.nl.error_message)
            self.ncp.add("nonlinear")
            timer.end("nonlinear")

        if "transfer" in level:
            timer.start("transfer")
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.nl), &(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)
            self.ncp.add("transfer")
            timer.end("transfer")

        if "spectra" in level:
            timer.start("spectra")
            if spectra_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.nl), &(self.tr),
                            &(self.sp)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sp.error_message)
            self.ncp.add("spectra")
            timer.end("spectra")

        if "lensing" in level:
            timer.start("lensing")
            if lensing_init(&(self.pr), &(self.pt), &(self.sp),
                            &(self.nl), &(self.le)) == _FAILURE_:
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

        timer.end("compute")

        if performance_report is not None:
            performance_report.update(timer.times)

        self.ready = True
        self.allocated = True

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
                Check and enforce the computation of the spectra module
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
        cdef double *rcl = <double*> calloc(self.sp.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md[index_md] = <double*> calloc(self.sp.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.sp.ct_size*self.sp.ic_ic_size[index_md], sizeof(double))

        # Define a list of integers, refering to the flags and indices of each
        # possible output Cl. It allows for a clear and concise way of looping
        # over them, checking if they are defined or not.
        has_flags = [
            (self.sp.has_tt, self.sp.index_ct_tt, 'tt'),
            (self.sp.has_ee, self.sp.index_ct_ee, 'ee'),
            (self.sp.has_te, self.sp.index_ct_te, 'te'),
            (self.sp.has_bb, self.sp.index_ct_bb, 'bb'),
            (self.sp.has_pp, self.sp.index_ct_pp, 'pp'),
            (self.sp.has_tp, self.sp.index_ct_tp, 'tp'),]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)

        if not spectra:
            raise CosmoSevereError("No Cl computed")
        lmaxR = self.sp.l_max_tot
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
            if spectra_cl_at_l(&self.sp, ell, rcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.sp.error_message)
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = rcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(rcl)
        for index_md in range(self.sp.md_size):
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
        cdef double *dcl = <double*> calloc(self.sp.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md[index_md] = <double*> calloc(self.sp.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.sp.ct_size*self.sp.ic_ic_size[index_md], sizeof(double))

        lmaxR = self.pt.l_lss_max
        has_flags = [
            (self.sp.has_dd, self.sp.index_ct_dd, 'dd'),
            (self.sp.has_td, self.sp.index_ct_td, 'td'),
            (self.sp.has_ll, self.sp.index_ct_ll, 'll'),
            (self.sp.has_dl, self.sp.index_ct_dl, 'dl'),
            (self.sp.has_tl, self.sp.index_ct_tl, 'tl')]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)
                l_max_flag = self.sp.l_max_ct[self.sp.index_md_scalars][index]
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
        size = (self.sp.d_size*(self.sp.d_size+1)-(self.sp.d_size-self.sp.non_diag)*
                (self.sp.d_size-1-self.sp.non_diag))/2;
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
            if spectra_cl_at_l(&self.sp, ell, dcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.sp.error_message)
            if 'dd' in spectra:
                for index in range(size):
                    cl['dd'][index][ell] = dcl[self.sp.index_ct_dd+index]
            if 'll' in spectra:
                for index in range(size):
                    cl['ll'][index][ell] = dcl[self.sp.index_ct_ll+index]
            if 'dl' in spectra:
                for index in range(size):
                    cl['dl'][index][ell] = dcl[self.sp.index_ct_dl+index]
            if 'td' in spectra:
                cl['td'][ell] = dcl[self.sp.index_ct_td]
            if 'tl' in spectra:
                cl['tl'][ell] = dcl[self.sp.index_ct_tl]
        cl['ell'] = np.arange(lmax+1)

        free(dcl)
        for index_md in range(self.sp.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def z_of_r (self,z_array):
        cdef double tau=0.0
        cdef int last_index=0 #junk
        cdef double * pvecback
        r = np.zeros(len(z_array),'float64')
        dzdr = np.zeros(len(z_array),'float64')

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        i = 0
        for redshift in z_array:
            if background_tau_of_z(&self.ba,redshift,&tau)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
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
        cdef double tau=0.0
        cdef int last_index = 0  # junk
        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba, z, &tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba, tau, self.ba.long_info,
                self.ba.inter_normal, &last_index, pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)
        lum_distance = pvecback[self.ba.index_bg_lum_distance]
        free(pvecback)
        return lum_distance

    # Gives the pk for a given (k,z)
    def pk(self,double k,double z):
        """
        Gives the pk for a given k and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check to verify if output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk
        cdef double pk_cb
        cdef double pk_velo
        cdef double pk_cross
        cdef int dummy

        # Quantities for the isocurvature modes
        cdef double *pk_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        cdef double *pk_cb_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. You must add mPk to the list of outputs."
                )

        if (self.nl.method == 0):
             if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,pk_ic,&pk_cb,pk_cb_ic)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)
        else:
             if spectra_pk_nl_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,&pk_cb) ==_FAILURE_:
                    raise CosmoSevereError(self.sp.error_message)

        free(pk_ic)
        free(pk_cb_ic)
        return pk

    # Gives the pk_cb for a given (k,z)
    def pk_cb(self,double k,double z):
        """
        Gives the pk_cb for a given k and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check to verify if output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk
        cdef double pk_cb
        cdef double pk_velo
        cdef double pk_cross
        cdef int dummy

        # Quantities for the isocurvature modes
        cdef double *pk_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        cdef double *pk_cb_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. You must add mPk to the list of outputs."
                )
        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use pk, rather than pk_cb."
                )

        if (self.nl.method == 0):
             if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,pk_ic,&pk_cb,pk_cb_ic)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)
        else:
             if spectra_pk_nl_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,&pk_cb) ==_FAILURE_:
                    raise CosmoSevereError(self.sp.error_message)

        free(pk_ic)
        free(pk_cb_ic)
        return pk_cb

    # Gives the linear pk for a given (k,z)
    def pk_lin(self,double k,double z):
        """
        Gives the linear pk for a given k and z (even if non linear corrections were requested to Class)

        .. note::

            there is an additional check to verify if output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk
        cdef double pk_cb
        cdef double pk_velo
        cdef double pk_cross
        cdef int dummy

        # Quantities for the isocurvature modes
        cdef double *pk_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        cdef double *pk_cb_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. You must add mPk to the list of outputs."
                )

        if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,pk_ic,&pk_cb,pk_cb_ic)==_FAILURE_:
            raise CosmoSevereError(self.sp.error_message)

        free(pk_ic)
        free(pk_cb_ic)
        return pk

    # Gives the linear pk for a given (k,z)
    def pk_cb_lin(self,double k,double z):
        """
        Gives the linear pk for a given k and z (even if non linear corrections were requested to Class)

        .. note::

            there is an additional check to verify if output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk
        cdef double pk_cb
        cdef double pk_velo
        cdef double pk_cross
        cdef int dummy

        # Quantities for the isocurvature modes
        cdef double *pk_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        cdef double *pk_cb_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. You must add mPk to the list of outputs."
                )
        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use pk_lin, rather than pk_cb_lin."
                )

        if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,pk_ic,&pk_cb,pk_cb_ic)==_FAILURE_:
            raise CosmoSevereError(self.sp.error_message)

        free(pk_ic)
        free(pk_cb_ic)
        return pk_cb

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

    # Gives sigma(R,z) for a given (R,z)
    def sigma(self,double R,double z):
        """
        Gives the pk for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs."
                )

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError(
                "In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc."
                )

        if spectra_sigma(&self.ba,&self.pm,&self.sp,R,z,&sigma)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)

        return sigma

    # Gives sigma_cb(R,z) for a given (R,z)
    def sigma_cb(self,double R,double z):
        """
        Gives the pk for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma_cb

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError(
                "No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs."
                )

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError(
                "In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc."
                )

        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use sigma, rather than sigma_cb."
                )

        if spectra_sigma_cb(&self.ba,&self.pm,&self.sp,R,z,&sigma_cb)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)

        return sigma_cb

    def age(self):
        self.compute(["background"])
        return self.ba.age

    def h(self):
        return self.ba.h

    def n_s(self):
        return self.pm.n_s

    def tau_reio(self):
        return self.th.tau_reio

    # Defined twice ?
    def Omega_m(self):
        return self.ba.Omega0_b+self.ba.Omega0_cdm+self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm

    # This is commented because in the current form it only applies
    # to minimal LambdaCDM.
    # On would need to add contributions from ncdm, ddmdr, etc.
    #def Omega_r(self):
    #    return self.ba.Omega0_g+self.ba.Omega0_ur

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

    def sigma8(self):
        self.compute(["spectra"])
        return self.sp.sigma8

    def sigma8_cb(self):
        self.compute(["spectra"])
        return self.sp.sigma8_cb

    def rs_drag(self):
        self.compute(["thermodynamics"])
        return self.th.rs_d

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return D_A

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D = pvecback[self.ba.index_bg_D]

        free(pvecback)

        return D

    def scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        f = pvecback[self.ba.index_bg_f]

        free(pvecback)

        return f

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        H = pvecback[self.ba.index_bg_H]

        free(pvecback)

        return H

    def ionization_fraction(self, z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,self.th.inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,self.th.inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
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

    def Omega0_m(self):
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        return self.ba.Omega0_b+self.ba.Omega0_cdm+self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm

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
        names = tmp.split("\t")[:-1]
        tmp = str(tmp.decode())
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


    @cython.returns(dict)
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.ccall
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
        conditions today. You must include 'dCl' and 'vCl' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that z is inside the region spanned by 'z_pk'.

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

        index_md = 0;
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if spectra_output_tk_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.op.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.sp.ln_k_size

        size_ic_data = timesteps*number_of_titles;
        ic_num = self.sp.ic_size[index_md];

        data = <double*>malloc(sizeof(double)*size_ic_data*ic_num)

        if spectra_output_tk_data(&self.ba, &self.pt, &self.sp, outf, <double> z, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.sp.error_message)

        spectra = {}

        for index_ic in range(ic_num):
            if spectra_firstline_and_ic_suffix(&self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
                raise CosmoSevereError(self.op.error_message)
            ic_key = <bytes> ic_suffix

            tmpdict = {}
            for i in range(number_of_titles):
                tmpdict[names[i]] = np.zeros(timesteps, dtype=np.double)
                for index in range(timesteps):
                    tmpdict[names[i]][index] = data[index_ic*size_ic_data+index*number_of_titles+i]

            if ic_num==1:
                spectra = tmpdict
            else:
                spectra[ic_key] = tmpdict

        free(titles)
        free(data)

        return spectra


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
                value = (self.ba.Omega0_b + self.ba.Omega0_cdm+
                         self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm)
            elif name == 'omega_m':
                value = (self.ba.Omega0_b + self.ba.Omega0_cdm+
                         self.ba.Omega0_ncdm_tot + self.ba.Omega0_dcdm)/self.ba.h**2
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
            # NEW FOR NEURAL NET
            elif name == 'rd_rec':
                value = self.th.rd_rec
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
            elif name == 'alpha_kp':
                value = self.sp.alpha_kp
            elif name == 'alpha_k1':
                value = self.sp.alpha_k1
            elif name == 'alpha_k2':
                value = self.sp.alpha_k2
            elif name == 'alpha_II_2_20':
                value = self.sp.alpha_II_2_20
            elif name == 'alpha_RI_2_20':
                value = self.sp.alpha_RI_2_20
            elif name == 'alpha_RR_2_20':
                value = self.sp.alpha_RR_2_20
            elif name == 'alpha_II_21_200':
                value = self.sp.alpha_II_21_200
            elif name == 'alpha_RI_21_200':
                value = self.sp.alpha_RI_21_200
            elif name == 'alpha_RR_21_200':
                value = self.sp.alpha_RR_21_200
            elif name == 'alpha_II_201_2500':
                value = self.sp.alpha_II_201_2500
            elif name == 'alpha_RI_201_2500':
                value = self.sp.alpha_RI_201_2500
            elif name == 'alpha_RR_201_2500':
                value = self.sp.alpha_RR_201_2500
            elif name == 'alpha_II_2_2500':
                value = self.sp.alpha_II_2_2500
            elif name == 'alpha_RI_2_2500':
                value = self.sp.alpha_RI_2_2500
            elif name == 'alpha_RR_2_2500':
                value = self.sp.alpha_RR_2_2500
            elif name == 'sigma8':
                value = self.sp.sigma8
            elif name == 'sigma8_cb':
                value = self.sp.sigma8_cb
            elif name == 'g_sd':
                value = self.sd.sd_parameter_table[0]
            elif name == 'y_sd':
                value = self.sd.sd_parameter_table[1]
            elif name == 'mu_sd':
                value = self.sd.sd_parameter_table[2]
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
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return k_nl

    def nonlinear_scale_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
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
        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use nonlinear_scale, rather than nonlinear_scale_cb."
                )
        for index_z in range(z_size):
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return k_nl_cb

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
        cdef int nonlinearint
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')
        nonlinearint=1 if nonlinear else 0
        spectra_fast_pk_at_kvec_and_zvec(&self.ba, &self.sp, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data, nonlinearint)
        return pk

    def get_pk_cb_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef int nonlinearint
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')
        nonlinearint=1 if nonlinear else 0
        spectra_fast_pk_at_kvec_and_zvec(&self.ba, &self.sp, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data, nonlinearint)
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
        corresponding conformal time values.
        """
        cdef:
            int index_tau;
            int bt_size = self.ba.bt_size;
            int bg_size = self.ba.bg_size;
            int index_bg_rs = self.ba.index_bg_rs;
            double * background_table = self.ba.background_table;
            double [:] r_s = np.zeros((bt_size));
            double [:] rho_b = np.zeros((bt_size));
            double [:] rho_g = np.zeros((bt_size));
            double [:] tau_bg = np.zeros((bt_size));
            double [:] a = np.zeros((bt_size));
            double [:] H = np.zeros((bt_size));
            double [:] D = np.zeros((bt_size));

        for index_tau in range(bt_size):
            tau_bg[index_tau] = self.ba.tau_table[index_tau]

            r_s[index_tau] = background_table[index_tau*bg_size+index_bg_rs]
            rho_b[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_rho_b];
            rho_g[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_rho_g];
            a[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_a];
            H[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_H];
            D[index_tau] = background_table[index_tau * bg_size + self.ba.index_bg_D];

        return {
                "tau": np.asarray(tau_bg),
                "r_s": np.asarray(r_s),
                "rho_b": np.asarray(rho_b),
                "rho_g": np.asarray(rho_g),
                "a": np.asarray(a),
                "H": np.asarray(H),
                "D": np.asarray(D),
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
            double [:] numpy_r_d = np.zeros((tt_size));
            double [:] numpy_g = np.zeros((tt_size));
            double [:] numpy_g_reco = np.zeros(tt_size);
            double [:] numpy_g_reio = np.zeros(tt_size);
            double [:] numpy_dg = np.zeros((tt_size));
            double [:] numpy_dg_reco = np.zeros((tt_size));
            double [:] numpy_dg_reio = np.zeros((tt_size));
            double [:] numpy_e_kappa = np.zeros((tt_size));
            double [:] numpy_dkappa = np.zeros((tt_size));
            double [:] numpy_tau = np.zeros((tt_size));

        for index_z in range(tt_size):

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
            numpy_e_kappa[index_z] = thermodynamics_table[index_z*th_size + index_th_exp_m_kappa]

        g_reco = np.asarray(numpy_g_reco)
        g_reio = np.asarray(numpy_g_reio)
        g_reco_prime = np.asarray(numpy_dg_reco)
        g_reio_prime = np.asarray(numpy_dg_reio)
        return {
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
                }

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D = pvecback[self.ba.index_bg_D]

        free(pvecback)

        return D

    def scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        f = pvecback[self.ba.index_bg_f]

        free(pvecback)

        return f

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        H = pvecback[self.ba.index_bg_H]

        free(pvecback)

        return H

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
        if background_at_tau(&self.ba,tau,self.ba.short_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
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

        # names = np.array([
        #     't0','t0_sw', 't0_isw',
        #     't1',
        #     't2','t2_reco', 't2_reio',
        #     'p',
        #     'delta_m','delta_g','theta_m','phi','phi_plus_psi','phi_prime','psi'
        #     ])

        names = []

        for index_k in range(k_size):
            k_array[index_k] = k[index_k]
        for index_tau in range(tau_size):
            tau_array[index_tau] = tau[index_tau]

        indices = []

        if self.pt.has_source_t:
            indices.extend([
                self.pt.index_tp_t0, self.pt.index_tp_t0_sw, self.pt.index_tp_t0_isw,
                self.pt.index_tp_t0_reco, self.pt.index_tp_t0_reio,
                self.pt.index_tp_t0_reco_no_isw, self.pt.index_tp_t0_reio_no_isw,
                self.pt.index_tp_t1,
                self.pt.index_tp_t2, self.pt.index_tp_t2_reco, self.pt.index_tp_t2_reio
                ])
            names.extend([
                "t0", "t0_sw", "t0_isw",
                "t0_reco", "t0_reio",
                "t0_reco_no_isw", "t0_reio_no_isw",
                "t1",
                "t2", "t2_reco", "t2_reio"
                ])
        if self.pt.has_source_p:
            indices.append(self.pt.index_tp_p)
            names.append("p")
        if self.pt.has_source_delta_m:
            indices.append(self.pt.index_tp_delta_m)
            names.append("delta_m")
        if self.pt.has_source_delta_g:
            indices.append(self.pt.index_tp_delta_g)
            names.append("delta_g")
        if self.pt.has_source_theta_m:
            indices.append(self.pt.index_tp_theta_m)
            names.append("theta_m")
        if self.pt.has_source_phi:
            indices.append(self.pt.index_tp_phi)
            names.append("phi")
        if self.pt.has_source_phi_plus_psi:
            indices.append(self.pt.index_tp_phi_plus_psi)
            names.append("phi_plus_psi")
        if self.pt.has_source_phi_prime:
            indices.append(self.pt.index_tp_phi_prime)
            names.append("phi_prime")
        if self.pt.has_source_psi:
            indices.append(self.pt.index_tp_psi)
            names.append("psi")

        for index_type, name in zip(indices, names):
            for index_k in range(k_size):
                for index_tau in range(tau_size):
                    tmparray[index_k][index_tau] = sources_ptr[index_md][index_ic*tp_size+index_type][index_tau*k_size + index_k];

            sources[name] = np.asarray(tmparray)
            tmparray = np.zeros((k_size,tau_size))

        return (sources, np.asarray(k_array), np.asarray(tau_array))


