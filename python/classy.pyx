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

# TODO conditional import?
# TODO relative import?
import classynet.workspace
import classynet.predictors

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

cdef class c_linked_list:
    cdef clist_node* tail

    def __cinit__(self):
      # The list is empty
      self.tail = NULL

    def __dealloc__(self):
      while not (self.tail==NULL):
        temp = self.tail.prev
        free(self.tail)
        self.tail=temp

    cdef append(self,item):
      # Allocate new node
      newnode = <clist_node*>malloc(sizeof(clist_node))
      newnode.next = NULL
      newnode.prev = self.tail
      # Setting the value is a bit difficult, and proceeds in two steps:
      sanitized_item = string_to_char(item) # First, make sure it's a bytes object
      strcpy(newnode.value,sanitized_item)  # Second, copy the string into the node (at most lenght 40)
      # Connect previous node with new node, if that exists
      if self.tail != NULL:
        self.tail.next = newnode
      # Update what is the 'tail node'=last node
      self.tail = newnode

    cdef pop(self):
      # Keep temporary reference to 'tail node'=last node
      temp = self.tail
      # Retrieve the memory stored in it
      content = temp.value
      # Now we can discard the node
      self.tail = temp.prev
      if temp.prev != NULL:
        temp.prev.next = NULL
      free(temp)
      # Now we can return the memory
      return content

    cdef empty(self):
      return (self.tail==NULL)

    # Check if the linked list contains a single value
    cdef contains(self,value):
      temp = self.tail
      sanitized_value = string_to_char(value)
      while temp!=NULL:
        if(temp.value==sanitized_value):
          return True
        temp = temp.prev
      return False

    # Check if the linked list contains all values
    cdef contains_all(self,values):
      flag = True
      for value in values:
        if not self.contains(value):
          flag=False
          break
      return flag

    cdef clean(self):
      while not self.empty():
        self.pop()


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

    cpdef int computed # Flag to see if classy has already computed with the given pars
    cpdef int allocated # Flag to see if classy structs are allocated already
    cpdef object _pars # Dictionary of the parameters
    cpdef object ncp   # Keeps track of the structures initialized, in view of cleaning.
    cdef c_linked_list module_list

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
            return True
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
        self.allocated = False
        self.use_NN = False
        self.computed = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.module_list = c_linked_list()
        if default: self.set_default()

    def __dealloc__(self):
        if self.allocated:
          self.struct_cleanup()
        self.empty()
        # This part should always be done
        free(self.fc.filename)
        # Reset all the fc to zero if its not already done
        if self.fc.size !=0:
            self.fc.size=0
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)

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

        # TODO if pars contains NN enable, load models here?
        # TODO if pars contains NN enable, set use_NN

        return True

    def use_nn(self):
        """
        Utility methods that returns whether neural networks are enabled
        by checking whether 'neural network path' is in the input parameters.
        Also checks the 'nn_verbose' parameter to determine how much information 
        to print about the usage of neural networks.
        """
        if not "neural network path" in self._pars:
            return False
        elif "neural network path" in self.pars:
            using_nn = self.can_use_nn()
            if "nn_verbose" in self.pars:
                if self.pars["nn_verbose"]>1:
                    if not using_nn:
                        print("##################################")
                        print("#   NOT USING NEURAL NETWORKS!   #")
                        print("##################################")
                        return False
                    else:
                        print("##################################")
                        print("#    USING NEURAL NETWORKS!      #")
                        print("##################################")
                        return True
                elif self.pars["nn_verbose"]==1:
                    if not using_nn:
                        print("USING NEURAL NETWORKS : False!")
                        return False
                    else:
                        print("USING NEURAL NETWORKS :          True!")
                        return True
                elif self.pars["nn_verbose"]==0:
                    if not using_nn:
                        return False
                    else:
                        return True
                else:
                    raise ValueError("nn_verbose is not set to valid value: should be integer above 0 but is {}".format(self.pars["nn_verbose"]))
        


    def can_use_nn(self):
        """ may only be called if neural networks are enabled """
        workspace = self.nn_workspace()
        domain = workspace.loader().domain_descriptor()
        using_nn, self.pt.network_deltachisquared = domain.contains(self._pars)
        if not using_nn:
            if "nn_verbose" in self.pars:
                if self.pars["nn_verbose"]>1:
                    print("neural network domain of validity does not contain requested parameters")
            else:
                print("neural network domain of validity does not contain requested parameters")
            return False

        def expect(key, value):
            if not key in self._pars:
                print("expected key '{}' not found in parameters.".format(key))
                return False
            else:
                found = self._pars[key]
                if found != value:
                    print("expected parameter '{}' to be {}; got {} instead.".format(key, value, found))
                    return False
                else:
                    return True

        if not expect("N_ncdm", 1):
            return False
        if not expect("deg_ncdm", 3):
            return False
        if not expect("Omega_Lambda", 0):
            return False
        # TODO are there other valid values (e.g. 'true' or something like that)?
        if not expect("compute damping scale", "yes"):
            return False

        pk_max = self._pars.get("P_k_max_1/Mpc")
        if pk_max is not None and pk_max > 100.0:
            print("neural networks only applicable with 'P_k_max_1/Mpc' <= 100.0")
            return False

        return True

    def empty(self):
        self._pars = {}
        self.computed = False

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
        if self.module_list.contains("distortions"):
            distortions_free(&self.sd)
        if self.module_list.contains("lensing"):
            lensing_free(&self.le)
        if self.module_list.contains("spectra"):
            spectra_free(&self.sp)
        if self.module_list.contains("transfer"):
            transfer_free(&self.tr)
        if self.module_list.contains("nonlinear"):
            nonlinear_free(&self.nl)
        if self.module_list.contains("primordial"):
            primordial_free(&self.pm)
        if self.module_list.contains("perturb"):
            perturb_free(&self.pt)
        if self.module_list.contains("thermodynamics"):
            thermodynamics_free(&self.th)
        if self.module_list.contains("background"):
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
            if "spectra" not in level:
                level.append("spectra")
        if "spectra" in level:
            if "transfer" not in level:
                level.append("transfer")
        if "transfer" in level:
            if "nonlinear" not in level:
                level.append("nonlinear")
        if "nonlinear" in level:
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

    # GS: TODO remove this function
    def debug_overwrite_source(self, name, double[:, :] S):
        cdef int index_md = self.pt.index_md_scalars
        cdef int index_ic = self.pt.index_ic_ad
        cdef int tp_size = self.pt.tp_size[index_md]
        cdef int tau_size = self.pt.tau_size
        cdef int k_size = self.pt.k_size[index_md]
        cdef int i_tau
        cdef int i_k

        index_type = self.translate_source_to_index(name)

        print("expected S.shape of", (k_size, tau_size))
        print("got      S.shape of", S.shape)

        for i_tau in range(tau_size):
            for i_k in range(k_size):
                self.pt.sources[index_md][index_ic*tp_size + index_type][i_tau*k_size + i_k] = S[i_k][i_tau]


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
        # if no other modules were requested, i.e. if self.module_list contains (or is
        # equivalent to) level. If it is the case, simply stop the execution of
        # the function.
        if self.computed and self.module_list.contains_all(level):
            return

        # Check if already allocated to prevent memory leaks
        if self.allocated:
            self.struct_cleanup()

        # Otherwise, proceed with the normal computation.
        self.computed = False

        # Equivalent of writing a parameter file
        self._fillparfile()

        # self.module_list will contain the list of computed modules
        self.module_list.clean()

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
                                    &self.pt, &self.tr, &self.pm, &self.sp,
                                    &self.nl, &self.le, &self.sd, &self.op, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.module_list.append("input")
            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            # GS: added this because neural network arguments are not relevant
            # to the C code.
            problematic_exceptions = set(["neural network path", "nn_cheat", "nn_debug","nn_verbose"])
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
            self.module_list.append("background")
            timer.end("background")

        if "thermodynamics" in level:
            timer.start("thermodynamics")
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.th.error_message)
            self.module_list.append("thermodynamics")
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
            double [:, :, :] NN_prediction
            double * c_NN_sources


        if "perturb" in level:

            timer.start("perturb")
            timer.start("perturb_init")


            # Allocate memory for ALL source functions (since transfer.c iterates over them)

            use_nn = self.use_nn()

            if use_nn and not self.nn_cheat_enabled():
                if "nn_verbose" in self.pars:
                    if self.pars["nn_verbose"]>2:
                        print("Using neural networks; skipping regular perturbation module.")
                self.pt.perform_NN_skip = _TRUE_

            if perturb_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            
            self.module_list.append("perturb")
            timer.end("perturb_init")

            # flag for using NN
            if use_nn:
                timer.start("neural network complete")
                timer.start("neural network initialization")

                index_md = self.pt.index_md_scalars;
                k_size = self.pt.k_size[index_md];
                tau_size = self.pt.tau_size;
                index_ic = self.pt.index_ic_ad;

                tp_size = self.pt.tp_size[index_md];

                tau_CLASS = np.zeros((tau_size))
                for index_tau in range(tau_size):
                    tau_CLASS[index_tau] = self.pt.tau_sampling[index_tau]

                requested_index_types = []

                # add all sources that class calculates to the list for predicting
                source_names = []
                if self.pt.has_source_t:
                    requested_index_types.extend([
                        self.pt.index_tp_t0,
                        self.pt.index_tp_t1,
                        self.pt.index_tp_t2,
                        ])
                    source_names.extend(["t0", "t1", "t2"])

                    if self.nn_debug_enabled():
                        requested_index_types.extend([
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
                if self.pt.has_source_phi_plus_psi:
                    requested_index_types.append(self.pt.index_tp_phi_plus_psi)
                    source_names.append("phi_plus_psi")

                if self.pt.has_source_delta_m:
                    requested_index_types.append(self.pt.index_tp_delta_m)
                    source_names.append("delta_m")


                if self.pt.has_source_delta_cb:
                    requested_index_types.append(self.pt.index_tp_delta_cb)
                    source_names.append("delta_cb")

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
                timer.start("build predictor")
                predictor = classynet.predictors.build_predictor(self)
                timer.end("build predictor")
                timer.start("predictor.predict_many")
                k_NN, NN_prediction = predictor.predict_many(source_names, np.asarray(tau_CLASS))
                timer.end("predictor.predict_many")
                timer.end("get all sources")

                timer.start("overwrite k array")

                # Copy k values from NN
                k_NN_size = len(k_NN)
                free(self.pt.k[index_md])
                self.pt.k[index_md] = <double*>malloc(k_NN_size * sizeof(double))
                for i_k in range(k_NN_size):
                    self.pt.k[index_md][i_k] = k_NN[i_k]

                self.pt.k_min = k_NN[0]
                self.pt.k_max = k_NN[-1]
                self.pt.k_size[index_md] = k_NN_size

                k_max_cl = 0.4
                # TODO - 1 maybe?
                k_max_cl_idx = np.searchsorted(k_NN, k_max_cl)
                # self.pt.k_size_cl[index_md] = k_NN_size
                self.pt.k_size_cl[index_md] = k_max_cl_idx

                _k_max_dbg = self.pt.k[index_md][self.pt.k_size_cl[index_md] - 1]
                if "nn_verbose" in self.pars:
                    if self.pars["nn_verbose"]>2:
                        print("pt.k[index_md][pt.k_size_cl[index_md] - 1] =", _k_max_dbg)


                timer.end("overwrite k array")

                timer.start("allocate unused source functions")

                for index_type in range(tp_size):
                    # Using malloc instead of calloc here will cause the splining
                    # in transfer.c to explode, but that doesn't seem to be an issue.
                    # Using malloc over calloc saves about a factor of 10 in runtime.
                    # self.pt.sources[index_md][index_ic * tp_size + index_type] = <double*> calloc(k_NN_size * tau_size,  sizeof(double))
                    self.pt.sources[index_md][index_ic * tp_size + index_type] = <double*> malloc(k_NN_size * tau_size * sizeof(double))
                timer.end("allocate unused source functions")

                for key, value in predictor.times.items():
                    timer[key] = value

                for key, value in predictor.time_prediction_per_network.items():
                    timer["indiv. network: '{}'".format(key)] = value

                timer.start("overwrite source functions")

                for i, index_tp_x in enumerate(requested_index_types):
                    self.overwrite_source_function(
                            index_md, index_ic,
                            index_tp_x,
                            k_NN_size, tau_size, NN_prediction[i, :, :]
                            )
                #if self.pt.has_source_delta_cb:
                #    self.overwrite_source_function(
                #            index_md, index_ic,
                #            self.pt.index_tp_delta_cb,
                #            k_NN_size, tau_size, NN_prediction[source_names.index("delta_m"), :, :]
                #            )


                ############################################################
                # if self.pt.has_source_t:
                #     self.overwrite_source_function(
                #             index_md, index_ic,
                #             self.pt.index_tp_t0,
                #             k_NN_size, tau_size, NN_prediction[0, :, :]
                #             )
                #     self.overwrite_source_function(
                #             index_md, index_ic,
                #             self.pt.index_tp_t1,
                #             k_NN_size, tau_size, NN_prediction[1, :, :]
                #             )
                #     self.overwrite_source_function(
                #             index_md, index_ic,
                #             self.pt.index_tp_t2,
                #             k_NN_size, tau_size, NN_prediction[2, :, :]
                #             )

                if self.pt.has_source_p:
                    assert "t2" in source_names
                    self.overwrite_source_function(
                        index_md, index_ic,
                        self.pt.index_tp_p,
                        k_NN_size, tau_size,
                        np.sqrt(6) * NN_prediction[source_names.index("t2"), :, :]
                    )

                # if self.pt.has_source_phi_plus_psi:
                #     self.overwrite_source_function(
                #             index_md, index_ic,
                #             self.pt.index_tp_phi_plus_psi,
                #             k_NN_size, tau_size, NN_prediction[3, :, :]
                #             )

                # if self.pt.has_source_delta_m:
                #     self.overwrite_source_function(
                #             index_md, index_ic,
                #             self.pt.index_tp_delta_m,
                #             k_NN_size, tau_size, NN_prediction[4, :, :]
                #             )


                timer.end("overwrite source functions")
                ############################################################

                timer.end("neural network complete")

            timer.end("perturb")

            #print(self.ncp)

            if post_perturb_callback:
                post_perturb_callback(self)


        if "primordial" in level:
            timer.start("primordial")
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            
            self.module_list.append("primordial")
            timer.end("primordial")


        if "nonlinear" in level:
            timer.start("nonlinear")
            if nonlinear_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.nl) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.nl.error_message)

            self.module_list.append("nonlinear")
            timer.end("nonlinear")


        if "transfer" in level:
            timer.start("transfer")
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.nl), &(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)

            self.module_list.append("transfer")
            timer.end("transfer")

        if "spectra" in level:
            timer.start("spectra")
            if spectra_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.nl), &(self.tr),
                            &(self.sp)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sp.error_message)
            self.module_list.append("spectra")
            timer.end("spectra")


        if "lensing" in level:
            timer.start("lensing")
            if lensing_init(&(self.pr), &(self.pt), &(self.sp),
                            &(self.nl), &(self.le)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.le.error_message)

            self.module_list.append("lensing")
            timer.end("lensing")

        if "distortions" in level:
            timer.start("distortions")
            if distortions_init(&(self.pr), &(self.ba), &(self.th),
                                &(self.pt), &(self.pm), &(self.sd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sd.error_message)
            
            self.module_list.append("distortions")
            timer.end("distortions")


        timer.end("compute")

        if performance_report is not None:
            performance_report.update(timer.times)

        self.computed = True
        #print(self.ncp)
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
        size = int((self.sp.d_size*(self.sp.d_size+1)-(self.sp.d_size-self.sp.non_diag)*
                (self.sp.d_size-1-self.sp.non_diag))/2);
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

            if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
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

        if background_at_tau(&self.ba, tau, long_info,
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

        if (self.nl.method == nl_none):
            #print("no nonlinear_method called_here")
            if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_m,&pk,NULL)==_FAILURE_:
                raise CosmoSevereError(self.nl.error_message)
        else:
            if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_nonlinear,k,z,self.nl.index_pk_m,&pk,NULL)==_FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return pk

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
        if (self.nl.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        if (self.nl.method == nl_none):
            if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.nl.error_message)
        else:
            if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_nonlinear,k,z,self.nl.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

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

        if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_m,&pk_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

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

        if (self.nl.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_cb,&pk_cb_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

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


    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species = False):
        """
        Returns a grid of matter power spectrum values and the z and k
        at which it has been fully computed. Useful for creating interpolators.

        Parameters
        ----------
        nonlinear : bool
                Whether the returned power spectrum values are linear or non-linear (default)
        nonlinear : bool
                Whether the returned power spectrum is for galaxy clustering and excludes massive neutrinos, or always includes evrything (default)
        """
        cdef np.ndarray[DTYPE_t,ndim=2] pk_at_k_z = np.zeros((self.nl.k_size, self.nl.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.nl.k_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.nl.ln_tau_size),'float64')
        cdef int index_k, index_tau, index_pk
        cdef double z_max_nonlinear, z_max_requested

        # consistency checks

        if self.nl.has_pk_matter == False:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations; add 'mPk' in 'output'")

        if nonlinear == True and self.nl.method == nl_none:
            raise CosmoSevereError("You ask classy to return an array of nonlinear P(k,z) values, but the input parameters sent to CLASS did not require any non-linear P(k,z) calculations; add e.g. 'halofit' or 'HMcode' in 'nonlinear'")

        # check wich type of P(k) to return (total or clustering only, i.e. without massive neutrino contribution)
        if (only_clustering_species == True):
            index_pk = self.nl.index_pk_cluster
        else:
            index_pk = self.nl.index_pk_total

        # get list of redshfits

        if self.nl.ln_tau_size == 1:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in xrange(self.nl.ln_tau_size):
                if index_tau == self.nl.ln_tau_size-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.nl.ln_tau[index_tau]))

        # check consitency of the list of redshifts

        if nonlinear == True:
            z_max_nonlinear = self.z_of_tau(self.nl.tau[self.nl.index_tau_min_nl])
            z_max_requested = z[0]
            if ((self.nl.tau_size - self.nl.ln_tau_size) < self.nl.index_tau_min_nl):
                raise CosmoSevereError("get_pk_and_k_and_z() is trying to return P(k,z) up to z_max=%e (to encompass your requested maximum value of z); but the input parameters sent to CLASS were such that the non-linear P(k,z) could only be consistently computed up to z=%e; increase the input parameter 'P_k_max_h/Mpc' or 'P_k_max_1/Mpc', or increase the precision parameters 'halofit_min_k_max' and/or 'hmcode_min_k_max', or decrease your requested z_max"%(z_max_requested,z_max_nonlinear))

        # get list of k

        for index_k in xrange(self.nl.k_size):
            k[index_k] = self.nl.k[index_k]

        # get P(k,z) array

        for index_tau in xrange(self.nl.ln_tau_size):
            for index_k in xrange(self.nl.k_size):
                if nonlinear == True:
                    pk_at_k_z[index_k, index_tau] = np.exp(self.nl.ln_pk_nl[index_pk][index_tau * self.nl.k_size + index_k])
                else:
                    pk_at_k_z[index_k, index_tau] = np.exp(self.nl.ln_pk_l[index_pk][index_tau * self.nl.k_size + index_k])

        return pk_at_k_z, k, z

    # Gives sigma(R,z) for a given (R,z)
    def sigma(self,double R,double z):
        """
        Gives sigma (total matter) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        if nonlinear_sigmas_at_z(&self.pr,&self.ba,&self.nl,R,z,self.nl.index_pk_m,out_sigma,&sigma)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

        return sigma

    # Gives sigma_cb(R,z) for a given (R,z)
    def sigma_cb(self,double R,double z):
        """
        Gives sigma (cdm+b) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma_cb

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.nl.has_pk_cb == _FALSE_):
            raise CosmoSevereError("sigma_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        if nonlinear_sigmas_at_z(&self.pr,&self.ba,&self.nl,R,z,self.nl.index_pk_cb,out_sigma,&sigma_cb)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

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

        if (k < self.nl.k[1] or k > self.nl.k[self.nl.k_size-2]):
            raise CosmoSevereError("In order to get pk_tilt at k=%e 1/Mpc, you should compute P(k,z) in a wider range of k's"%k)

        if nonlinear_pk_tilt_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_total,&pk_tilt)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

        return pk_tilt

    #calculates the hmcode window_function of the Navarrow Frenk White Profile
    def nonlinear_hmcode_window_nfw(self,double k,double rv,double c):
        """
        Gives window_nfw for a given wavevector k, virial radius rv and concentration c

        """
        cdef double window_nfw


        if nonlinear_hmcode_window_nfw(&self.nl,k,rv,c,&window_nfw)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)

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
        self.compute(["nonlinear"])
        return self.nl.sigma8[self.nl.index_pk_m]

    #def neff(self):
    #    self.compute(["spectra"])
    #    return self.sp.neff

    def sigma8_cb(self):
        self.compute(["nonlinear"])
        return self.nl.sigma8[self.nl.index_pk_cb]

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return D_A

    def comoving_distance(self, z):
        """
        comoving_distance(z)

        Return the comoving distance

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

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        r = pvecback[self.ba.index_bg_conf_distance]

        free(pvecback)

        return r

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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        Om_m = pvecback[self.ba.index_bg_Omega_m]

        free(pvecback)

        return Om_m


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

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
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
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
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

        if perturb_output_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pt.k_size[index_md]

        size_ic_data = timesteps*number_of_titles;
        ic_num = self.pt.ic_size[index_md];

        data = <double*>malloc(sizeof(double)*size_ic_data*ic_num)

        if perturb_output_data(&self.ba, &self.pt, outf, <double> z, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        transfers = {}

        for index_ic in range(ic_num):
            if perturb_output_firstline_and_ic_suffix(&self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
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
                value = self.nl.sigma8[self.nl.index_pk_m]
            elif name == 'sigma8_cb':
                value = self.nl.sigma8[self.nl.index_pk_cb]
            elif name == 'k_eq':
                value = self.ba.a_eq*self.ba.H_eq
            elif name == 'g_sd':
                value = self.sd.sd_parameter_table[0]
            elif name == 'y_sd':
                value = self.sd.sd_parameter_table[1]
            elif name == 'mu_sd':
                value = self.sd.sd_parameter_table[2]
            elif name == 'network_delta_chi_sqared':
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
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

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
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return k_nl_cb

    def nonlinear_hmcode_sigma8(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigma8_at_z(&self.ba,&self.nl,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_8

    def nonlinear_hmcode_sigma8_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigma8_at_z(&self.ba,&self.nl,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_8_cb

    def nonlinear_hmcode_sigmadisp(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp(z, z_size)

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
#            if nonlinear_hmcode_sigmadisp_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp

    def nonlinear_hmcode_sigmadisp_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp(z, z_size)

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
#            if nonlinear_hmcode_sigmadisp_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp_cb

    def nonlinear_hmcode_sigmadisp100(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp100(z, z_size)

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
#            if nonlinear_hmcode_sigmadisp100_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp_100

    def nonlinear_hmcode_sigmadisp100_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp100(z, z_size)

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
#            if nonlinear_hmcode_sigmadisp100_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp_100_cb

    def nonlinear_hmcode_sigmaprime(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmaprime(z, z_size)

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
#            if nonlinear_hmcode_sigmaprime_at_z(&self.ba,&self.nl,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_prime

    def nonlinear_hmcode_sigmaprime_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmaprime(z, z_size)

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
#            if nonlinear_hmcode_sigmaprime_at_z(&self.ba,&self.nl,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

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
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk

    def get_pk_cb_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

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
            numpy_e_kappa[index_z] = thermodynamics_table[index_z*th_size + index_th_exp_m_kappa]

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
                }
        return ret

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

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
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

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
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

        if background_at_z(&self.ba,z,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
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
        workspace = self._pars["neural network path"]
        if any(isinstance(workspace, t) for t in [str, bytes, os.PathLike]):
            workspace = classynet.workspace.Workspace(workspace)
        return workspace

    def nn_cosmological_parameters(self):
        manifest = self.nn_workspace().loader().manifest()
        names = manifest["cosmological_parameters"]

        result = {}
        remaining = []
        for name in names:
            if name in self._pars:
                result[name] = self._pars[name]
            else:
                remaining.append(name)

        for name in remaining:
            if name == "omega_b":
                result[name] = self.omega_b()
            elif name == "omega_cdm":
                result[name] = self.omega_cdm()
            elif name == "h":
                result[name] = self.h()
            elif name == "tau_reio":
                result[name] = self.tau_reio()
            elif name == "Omega_k":
                result[name] = self.get_current_derived_parameters(["Omega_k"])["Omega_k"]
            # Regarding w0_fld and wa_fld: It is verified that Omega_Lambda=0 in `can_use_nn`.
            elif name == "w0_fld":
                result[name] = 0.0
            elif name == "w0_fld":
                result[name] = -1.0
            else:
                raise ValueError("Unknown parameter: '{}'".format(name))

        return result

    def nn_cheat_enabled(self):
        return "nn_cheat" in self._pars

    def nn_cheat_sources(self):
        return self._pars["nn_cheat"]

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


    def _debug_transfer(self):
        k, tau = self.get_k_tau()
        kmintau0 = k[0] * tau[-1]
        print("tau0 = tau[-1] = {}".format(tau[-1]))
        print("k_min = k[0] = {}".format(k[0]))
        print("kmintau0 = k[0] * tau[-1] = {}".format(kmintau0))

        if self.ba.sgnK == 0:
            print("K = 0, flat case")
            # K<0 (flat)  : start close to zero
            print("kmintau0 / conformal_age =", kmintau0 / self.ba.conformal_age)
        elif self.ba.sgnK == -1:
            print("K = {} < 0, open case".format(self.ba.K))
            # K<0 (open)  : start close to sqrt(-K)
            # (in transfer modules, for scalars, this will correspond to q close to zero;
            # for vectors and tensors, this value is even smaller than the minimum necessary value)
            # return np.sqrt(-self.ba.K + pow(self.pr.k_min_tau0 / self.ba.conformal_age / self.th.angular_rescaling, 2))
        elif self.ba.sgnK == 1:
            print("K = {} > 0, closed case".format(self.ba.K))
            # K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K)
            # return np.sqrt((8.-1.e-4) * self.ba.K);

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
