# distutils: language = c++
"""
.. module:: classy
    :synopsis: Python wrapper around CLASS
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>
.. moduleauthor:: Thomas Tram <thomas.tram@phys.au.dk>

This module defines a class called Class. It is used with Monte Python to
extract cosmological parameters.

"""
DEF _FALSE_ = 0
DEF _FAILURE_ = 1
DEF _MAXTITLESTRINGLENGTH_ = 8000

from libc.math cimport exp, log
from libc.stdlib cimport malloc, calloc, free
from libc.stdio cimport sprintf

from libcpp cimport bool
from libcpp.memory cimport nullptr
from libcpp.memory cimport unique_ptr
from libcpp.memory cimport shared_ptr
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

from numpy.math cimport EULER, LOGE2

import cython
cimport cython
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np
import sys
from cclassy cimport *
# Nils : Added for python 3.x and python 2.x compatibility
cdef viewdictitems(dict d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=None):
        if message is None:
            message = ''
        elif isinstance(message, bytes):
            try:
                # If message has garbage inside, decode may fail.
                message = message.decode()
            except:
                pass

        self.message = message

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

cdef int raise_my_py_error() except *:
    cdef:
        pair[string, string] cpp_exception
        str cpp_exception_type
    cpp_exception = get_my_py_error_message()
    cpp_exception_type = cpp_exception.first.decode()
    if "invalid_argument" in cpp_exception_type:
        raise CosmoSevereError(cpp_exception.second)
    elif "runtime_error" in cpp_exception_type:
        raise CosmoComputationError(cpp_exception.second)
    else:
        raise NotImplementedError(cpp_exception.second)

cdef extern from "cosmology.h":
    ctypedef shared_ptr[const InputModule] InputModulePtr
    ctypedef shared_ptr[const BackgroundModule] BackgroundModulePtr
    ctypedef shared_ptr[const ThermodynamicsModule] ThermodynamicsModulePtr
    ctypedef shared_ptr[const PerturbationsModule] PerturbationsModulePtr
    ctypedef shared_ptr[const PrimordialModule] PrimordialModulePtr
    ctypedef shared_ptr[const NonlinearModule] NonlinearModulePtr
    ctypedef shared_ptr[const TransferModule] TransferModulePtr
    ctypedef shared_ptr[const SpectraModule] SpectraModulePtr
    ctypedef shared_ptr[const LensingModule] LensingModulePtr

    cdef cppclass Cosmology:
        Cosmology(FileContent& fc) except +raise_my_py_error
        InputModulePtr& GetInputModule()
        BackgroundModulePtr& GetBackgroundModule() except +raise_my_py_error
        ThermodynamicsModulePtr& GetThermodynamicsModule() except +raise_my_py_error
        PerturbationsModulePtr& GetPerturbationsModule() except +raise_my_py_error
        PrimordialModulePtr& GetPrimordialModule() except +raise_my_py_error
        NonlinearModulePtr& GetNonlinearModule() except +raise_my_py_error
        TransferModulePtr& GetTransferModule() except +raise_my_py_error
        SpectraModulePtr& GetSpectraModule() except +raise_my_py_error
        LensingModulePtr& GetLensingModule() except +raise_my_py_error

# To support the legacy name Class for the cosmology class.
cdef class Class(PyCosmology):
    pass

cdef class PyCosmology:

    """
    Cython wrapper class for C++ class Cosmology
    """

    cdef unique_ptr[Cosmology] _thisptr
    cdef dict _pars
    cdef bool parameters_changed
    cdef FileContent _fc

    cdef const precision* pr
    cdef const background* ba
    cdef const thermo* th
    cdef const perturbs* pt
    cdef const primordial* pm
    cdef const nonlinear* nl
    cdef const transfers* tr
    cdef const spectra* sp
    cdef const lensing* le
    cdef const output* op

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

    def __init__(self, input_parameters=None):
        if input_parameters is None:
            input_parameters = {}
        self._pars = input_parameters
        self.reset()

    cdef reset(self):
        cdef:
            Py_ssize_t i
            bool problem_flag
        self._update_fc_from_pars()
        self.parameters_changed = False
        self._thisptr.reset(new Cosmology(self._fc))
        # This part is done to list all the unread parameters, for debugging
        problem_flag = False
        problematic_parameters = []
        for i in range(self._fc.size):
            if self._fc.read[i] == _FALSE_:
                problem_flag = True
                problematic_parameters.append(self._fc.name[i].decode())
        if problem_flag:
            raise CosmoSevereError(
                "Class did not read input parameter(s): {}\n"
                .format(', '.join(problematic_parameters))
            )

        input_module = deref(self._thisptr).GetInputModule()
        self.pr = &deref(input_module).precision_
        self.ba = &deref(input_module).background_
        self.th = &deref(input_module).thermodynamics_
        self.pt = &deref(input_module).perturbations_
        self.pm = &deref(input_module).primordial_
        self.nl = &deref(input_module).nonlinear_
        self.tr = &deref(input_module).transfers_
        self.sp = &deref(input_module).spectra_
        self.le = &deref(input_module).lensing_
        self.op = &deref(input_module).output_
        return self

    cdef _update_fc_from_pars(self):
        cdef:
            FileContent new_file_content
            char* dumc
            int i

        # _fc will be cleaned up by its destructor
        self._fc = new_file_content
        self._fc.size = len(self._pars)
        self._fc.name = <FileArg*> malloc(sizeof(FileArg)*self._fc.size)
        assert(self._fc.name != nullptr)
        self._fc.value = <FileArg*> malloc(sizeof(FileArg)*self._fc.size)
        assert(self._fc.value != nullptr)
        self._fc.read = <short*> malloc(sizeof(short)*self._fc.size)
        assert(self._fc.read != nullptr)

        # fill parameter file
        for i, (name, value) in enumerate(self._pars.items()):
            dumcp = name.encode()
            dumc = dumcp
            sprintf(self._fc.name[i], "%s", dumc)

            dumcp = str(value).encode()
            dumc = dumcp
            sprintf(self._fc.value[i], "%s", dumc)
            self._fc.read[i] = 0
        return self

    # The functions struct_cleanup(), empty(), set() and compute() are not neccessary, but they are here to support
    # legacy code and MontePython.
    cpdef struct_cleanup(self):
        return

    cpdef empty(self):
        self._pars = {}
        self.parameters_changed = True
        return self

    cpdef set(self, input_parameters):
        if viewdictitems(input_parameters) <= viewdictitems(self._pars):
            return
        self._pars.update(input_parameters)
        self.parameters_changed = True
        return self

    cpdef compute(self, level=None):
        if level is None:
            level = ['lensing']
        if self.parameters_changed:
            self.reset()
        final_level = level[0].lower()
        if final_level == 'background':
            deref(self._thisptr).GetBackgroundModule()
        elif final_level == 'thermodynamics':
            deref(self._thisptr).GetThermodynamicsModule()
        elif final_level == 'perturb':
            deref(self._thisptr).GetPerturbationsModule()
        elif final_level == 'primordial':
            deref(self._thisptr).GetPrimordialModule()
        elif final_level == 'nonlinear':
            deref(self._thisptr).GetNonlinearModule()
        elif final_level == 'transfer':
            deref(self._thisptr).GetTransferModule()
        elif final_level == 'spectra':
            deref(self._thisptr).GetSpectraModule()
        elif final_level == 'lensing':
            deref(self._thisptr).GetLensingModule()
        return self

    cpdef get_input_precision(self):
        return deref(self.pr)

    cpdef get_input_background(self):
        return deref(self.ba)

    cpdef get_input_thermodynamics(self):
        return deref(self.th)

    cpdef get_input_perturbations(self):
        return deref(self.pt)

    cpdef get_input_transfers(self):
        return deref(self.tr)

    cpdef get_input_primordial(self):
        return deref(self.pm)

    cpdef get_input_spectra(self):
        return deref(self.sp)

    cpdef get_input_nonlinear(self):
        return deref(self.nl)

    cpdef get_input_lensing(self):
        return deref(self.le)

    cpdef get_input_output(self):
        return deref(self.op)

    cdef general_cl(self, int lmax, bool is_lensed):
        cdef:
            dict out_dict
            int lmax_computed
            int lmaxpp
            map[string, vector[double]] cl_data

        if self.pt.has_cls == _FALSE_:
            raise CosmoSevereError("No Cls computed")

        if is_lensed:
            le = deref(self._thisptr).GetLensingModule()
            lmax_computed = deref(le).l_lensed_max_
            if lmax == -1:
                lmax = lmax_computed
            if self.le.has_lensed_cls == _FALSE_:
                raise CosmoSevereError("Lensing Cls not computed, add 'lensing':'yes' to your input.")
            if lmax > lmax_computed:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmax_computed)
            cl_data = deref(le).cl_output(lmax)
        else:
            sp = deref(self._thisptr).GetSpectraModule()
            lmax_computed = deref(sp).l_max_tot_
            if lmax == -1:
                lmax = lmax_computed
            if lmax > lmax_computed:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmax_computed)
            cl_data = deref(sp).cl_output(lmax)
        lmaxpp = lmax + 1

        out_dict = {}
        for element in cl_data:
            key = <bytes> element.first
            key = str(key.decode())
            out_dict[key] = np.asarray(<double[:lmaxpp]> &element.second[0]).copy()
        out_dict['ell'] = np.arange(lmax + 1)
        return out_dict

    cpdef raw_cl_no_copy(self, int lmax=-1):
        cdef:
            dict out_dict
            double[::1] double_view
            int ct_size
            int lmax_computed
            map[string, int] index_map
            vector[double*] output_pointers

        sp = deref(self._thisptr).GetSpectraModule()
        lmax_computed = deref(sp).l_max_tot_
        if lmax == -1:
            lmax = lmax_computed
        if lmax > lmax_computed:
            raise CosmoSevereError("Can only compute up to lmax=%d"%lmax_computed)

        ct_size = deref(sp).ct_size_
        index_map = deref(sp).cl_output_index_map()
        output_pointers.resize(ct_size)
        out_dict = {}
        for element in index_map:
            key = <bytes> element.first
            key = str(key.decode())
            out_dict[key] = np.empty(lmax + 1, dtype=np.double)
            double_view = out_dict[key]
            output_pointers[element.second] = &double_view[0]

        deref(sp).cl_output_no_copy(lmax, output_pointers)
        out_dict['ell'] = np.arange(lmax + 1)
        return out_dict

    cpdef raw_cl(self, int lmax=-1):
        return self.general_cl(lmax, False)

    cpdef lensed_cl(self, int lmax=-1):
        return self.general_cl(lmax, True)

    cpdef z_of_r (self, double[::1] z_array):
        cdef:
            double tau
            int last_index
            int bg_size
            short long_info
            short inter_normal
            int index_bg_conf_distance
            int index_bg_H
            Py_ssize_t i
            Py_ssize_t z_array_size
            double z
            int status
            vector[double] pvecback
            double[::1] r
            double[::1] dzdr

        z_array_size = z_array.shape[0]
        r_arr = np.empty(z_array_size, np.double)
        r = r_arr
        dzdr_arr = np.empty(z_array_size, np.double)
        dzdr = dzdr_arr

        long_info = self.ba.long_info
        inter_normal = self.ba.inter_normal

        background_module = deref(self._thisptr).GetBackgroundModule()
        bg_size = deref(background_module).bg_size_
        index_bg_conf_distance = deref(background_module).index_bg_conf_distance_
        index_bg_H = deref(background_module).index_bg_H_

        pvecback.resize(bg_size)

        for i in range(z_array_size):
            z = z_array[i]
            status = deref(background_module).background_tau_of_z(z, &tau)
            if (status == _FAILURE_):
                raise CosmoSevereError(deref(background_module).error_message_)

            status = deref(background_module).background_at_tau(tau, long_info, inter_normal, &last_index, &pvecback[0])
            if (status == _FAILURE_):
                raise CosmoSevereError(deref(background_module).error_message_)

            # store r
            r[i] = pvecback[index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[index_bg_H]

        return r_arr, dzdr_arr

    cpdef luminosity_distance(self, double z):
        """
        luminosity_distance(z)
        """
        cdef int index_bg_lum_distance
        background_module = deref(self._thisptr).GetBackgroundModule()
        index_bg_lum_distance = deref(background_module).index_bg_lum_distance_
        return self.get_background_value_at_z(z, index_bg_lum_distance)

    cdef pk_general(self, double k, double z, int index_pk, pk_outputs linear_or_nonlinear):
        cdef:
            double pk
            int status

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("Power spectrum not computed. You must add mPk to the list of outputs.")

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        status = deref(nonlinear_module).nonlinear_pk_at_k_and_z(linear_or_nonlinear, k, z, index_pk, &pk, NULL)
        if status == _FAILURE_:
            raise CosmoSevereError(deref(nonlinear_module).error_message_)
        return pk

    cpdef pk(self, double k, double z):
        """
        Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef:
            int index_pk_m
            pk_outputs linear_or_nonlinear

        if self.nl.method == nl_none:
            linear_or_nonlinear = pk_linear
        else:
            linear_or_nonlinear = pk_nonlinear

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_m = deref(nonlinear_module).index_pk_m_

        return self.pk_general(k, z, index_pk_m, linear_or_nonlinear)

    # Gives the cdm+b pk for a given (k,z)
    cpdef pk_cb(self, double k, double z):
        """
        Gives the cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef:
            int index_pk_cb
            pk_outputs linear_or_nonlinear
            short has_pk_cb

        if self.nl.method == nl_none:
            linear_or_nonlinear = pk_linear
        else:
            linear_or_nonlinear = pk_nonlinear

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_cb = deref(nonlinear_module).index_pk_cb_
        has_pk_cb = deref(nonlinear_module).has_pk_cb_

        if (has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        return self.pk_general(k, z, index_pk_cb, linear_or_nonlinear)

    # Gives the total matter pk for a given (k,z)
    cpdef pk_lin(self, double k, double z):
        """
        Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef int index_pk_m
        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_m = deref(nonlinear_module).index_pk_m_

        return self.pk_general(k, z, index_pk_m, pk_linear)

    # Gives the cdm+b pk for a given (k,z)
    cpdef pk_cb_lin(self,double k,double z):
        """
        Gives the linear cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef:
            int index_pk_cb
            short has_pk_cb

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_cb = deref(nonlinear_module).index_pk_cb_
        has_pk_cb = deref(nonlinear_module).has_pk_cb_

        if (has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")
        return self.pk_general(k, z, index_pk_cb, pk_linear)

    cdef get_pk_general(self, double[:,:,::1] k, double[::1] z, Py_ssize_t k_size, Py_ssize_t z_size, Py_ssize_t mu_size, Py_ssize_t index_pk, pk_outputs linear_or_nonlinear):
        cdef:
            Py_ssize_t index_k
            Py_ssize_t index_z
            Py_ssize_t index_mu
            double pk_val
            double[:,:,::1] pk

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("Power spectrum not computed. You must add mPk to the list of outputs.")

        pk_arr = np.empty((k_size, z_size, mu_size), np.double)
        pk = pk_arr
        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        # TODO: Consider prange for outer loop
        for index_k in range(k_size):
            for index_z in range(z_size):
                for index_mu in range(mu_size):
                    pk[index_k, index_z, index_mu] = self.pk_general(k[index_k, index_z, index_mu], z[index_z], index_pk, linear_or_nonlinear)

        return pk_arr

    cpdef get_pk(self, double[:,:,::1] k, double[::1] z, Py_ssize_t k_size, Py_ssize_t z_size, Py_ssize_t mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef:
            int index_pk_m
            pk_outputs linear_or_nonlinear

        if self.nl.method == nl_none:
            linear_or_nonlinear = pk_linear
        else:
            linear_or_nonlinear = pk_nonlinear

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_m = deref(nonlinear_module).index_pk_m_
        return self.get_pk_general(k, z, k_size, z_size, mu_size, index_pk_m, linear_or_nonlinear)

    cpdef get_pk_cb(self, double[:,:,::1] k, double[::1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef:
            int index_pk_cb
            short has_pk_cb
            pk_outputs linear_or_nonlinear

        if self.nl.method == nl_none:
            linear_or_nonlinear = pk_linear
        else:
            linear_or_nonlinear = pk_nonlinear

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_cb = deref(nonlinear_module).index_pk_cb_
        has_pk_cb = deref(nonlinear_module).has_pk_cb_

        if (has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        return self.get_pk_general(k, z, k_size, z_size, mu_size, index_pk_cb, linear_or_nonlinear)

    cpdef get_pk_lin(self, double[:,:,::1] k, double[::1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef int index_pk_m
        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_m = deref(nonlinear_module).index_pk_m_

        return self.get_pk_general(k, z, k_size, z_size, mu_size, index_pk_m, pk_linear)

    cpdef get_pk_cb_lin(self, double[:,:,::1] k, double[::1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef:
            int index_pk_cb
            short has_pk_cb

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_cb = deref(nonlinear_module).index_pk_cb_
        has_pk_cb = deref(nonlinear_module).has_pk_cb_

        if (has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        return self.get_pk_general(k, z, k_size, z_size, mu_size, index_pk_cb, pk_linear)

    # Gives sigma(R,z) for a given (R,z)
    cpdef sigma(self, double R, double z):
        """
        Gives sigma (total matter) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef:
            double sigma
            int index_pk_m
            int status

        if self.pt.has_pk_matter == _FALSE_:
            raise CosmoSevereError("Power spectrum not computed. In order to get sigma(R, z) you must add mPk to the list of outputs.")
        if self.pt.k_max_for_pk < self.ba.h:
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_m = deref(nonlinear_module).index_pk_m_
        status = deref(nonlinear_module).nonlinear_sigmas_at_z(R, z, index_pk_m, out_sigma, &sigma)
        if status == _FAILURE_:
            raise CosmoSevereError(deref(nonlinear_module).error_message_)

        return sigma

    # Gives sigma(R, z) for a given (R, z)
    cpdef sigma_cb(self, double R, double z):
        """
        Gives sigma (cdm+b) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef:
            double sigma_cb
            int index_pk_cb
            int status

        if self.pt.has_pk_matter == _FALSE_:
            raise CosmoSevereError("Power spectrum not computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")
        if self.pt.k_max_for_pk < self.ba.h:
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        index_pk_cb = deref(nonlinear_module).index_pk_cb_
        has_pk_cb = deref(nonlinear_module).has_pk_cb_
        if (has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        status = deref(nonlinear_module).nonlinear_sigmas_at_z(R, z, index_pk_cb, out_sigma, &sigma_cb)
        if status == _FAILURE_:
            raise CosmoSevereError(deref(nonlinear_module).error_message_)

        return sigma_cb

    # Gives effective logarithmic slope of P_L(k,z) (total matter) for a given (k,z)
    cpdef pk_tilt(self, double k, double z):
        """
        Gives effective logarithmic slope of P_L(k,z) (total matter) for a given k and z
        (k is the wavenumber in units of 1/Mpc, z is the redshift, the output is dimensionless)

        .. note::

            there is an additional check to verify whether output contains `mPk` and whether k is in the right range

        """
        cdef:
            double ln_k
            double pk_tilt
            double* ln_k_vec
            int index_pk_m
            int k_size
            int status

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("Power spectrum not computed. In order to get pk_tilt(k, z) you must add mPk to the list of outputs.")

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        ln_k = log(k)
        ln_k_vec = deref(nonlinear_module).ln_k_
        k_size = deref(nonlinear_module).k_size_
        index_pk_m = deref(nonlinear_module).index_pk_m_
        if (k_size < 2 or not (ln_k_vec[1] <= ln_k <= ln_k_vec[k_size - 2])):
            raise CosmoSevereError("In order to get pk_tilt at k=%e 1/Mpc, you should compute P(k,z) in a wider range of k's"%k)
        status = deref(nonlinear_module).nonlinear_pk_tilt_at_k_and_z(pk_linear, k, z, index_pk_m, &pk_tilt)
        if status == _FAILURE_:
            raise CosmoSevereError(deref(nonlinear_module).error_message_)
        return pk_tilt


    cpdef age(self):
        background_module = deref(self._thisptr).GetBackgroundModule()
        return deref(background_module).age_

    cpdef h(self):
        return self.ba.h

    cpdef n_s(self):
        primordial_module = deref(self._thisptr).GetPrimordialModule()
        return deref(primordial_module).n_s_

    cpdef A_s(self):
        primordial_module = deref(self._thisptr).GetPrimordialModule()
        return deref(primordial_module).A_s_

    cpdef tau_reio(self):
        thm = deref(self._thisptr).GetThermodynamicsModule()
        return deref(thm).tau_reionization_

    cpdef Omega_m(self):
        background_module = deref(self._thisptr).GetBackgroundModule()
        return deref(background_module).Omega0_m_

    cpdef Omega_r(self):
        background_module = deref(self._thisptr).GetBackgroundModule()
        return deref(background_module).Omega0_r_

    cpdef theta_s_100(self):
        thm = deref(self._thisptr).GetThermodynamicsModule()
        return 100.*deref(thm).rs_rec_/deref(thm).ra_rec_

    cpdef theta_star_100(self):
        thm = deref(self._thisptr).GetThermodynamicsModule()
        return 100.*deref(thm).rs_star_/deref(thm).ra_star_

    cpdef Omega_Lambda(self):
        return self.ba.Omega0_lambda

    cpdef Omega_g(self):
        return self.ba.Omega0_g

    cpdef Omega_b(self):
        return self.ba.Omega0_b

    cpdef omega_b(self):
        return self.ba.Omega0_b*self.ba.h*self.ba.h

    cpdef Neff(self):
        bam = deref(self._thisptr).GetBackgroundModule()
        return deref(bam).Neff_

    cpdef k_eq(self):
        bam = deref(self._thisptr).GetBackgroundModule()
        return deref(bam).a_eq_*deref(bam).H_eq_

    cpdef sigma8(self):
        cdef int index_pk_m
        nlm = deref(self._thisptr).GetNonlinearModule()
        index_pk_m = deref(nlm).index_pk_m_
        return deref(nlm).sigma8_[index_pk_m]

    cpdef sigma8_cb(self):
        cdef int index_pk_cb
        nlm = deref(self._thisptr).GetNonlinearModule()
        index_pk_cb = deref(nlm).index_pk_cb_
        return deref(nlm).sigma8_[index_pk_cb]

    cpdef rs_drag(self):
        thm = deref(self._thisptr).GetThermodynamicsModule()
        return deref(thm).rs_d_

    cpdef z_reio(self):
        thm = deref(self._thisptr).GetThermodynamicsModule()
        return deref(thm).z_reionization_

    cdef get_background_value_at_z(self, double z, int index_bg):
        cdef:
            double output_value
            double tau
            int bg_size
            int last_index
            int status
            short inter_normal
            short long_info
            vector[double] pvecback

        long_info = self.ba.long_info
        inter_normal = self.ba.inter_normal

        background_module = deref(self._thisptr).GetBackgroundModule()
        bg_size = deref(background_module).bg_size_

        pvecback.resize(bg_size)
        status = deref(background_module).background_tau_of_z(z, &tau)
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)
        status = deref(background_module).background_at_tau(tau, long_info, inter_normal, &last_index, &pvecback[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)
        output_value = pvecback[index_bg]
        return output_value

    cdef get_thermodynamics_value_at_z(self, double z, int index_th):
        cdef:
            double tau
            int last_index
            vector[double] pvecback
            vector[double] pvecthermo
            int status
            short long_info
            short inter_normal
            int bg_size
            int th_size
            short inter_mode
            double output_value

        long_info = self.ba.long_info
        inter_normal = self.ba.inter_normal

        background_module = deref(self._thisptr).GetBackgroundModule()
        thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
        bg_size = deref(background_module).bg_size_
        th_size = deref(thermodynamics_module).th_size_

        pvecback.resize(bg_size)
        pvecthermo.resize(th_size)

        inter_mode = deref(thermodynamics_module).inter_normal_

        status = deref(background_module).background_tau_of_z(z, &tau)
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)
        status = deref(background_module).background_at_tau(tau, long_info, inter_normal, &last_index, &pvecback[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)

        status = deref(thermodynamics_module).thermodynamics_at_z(z, inter_mode, &last_index, &pvecback[0], &pvecthermo[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(thermodynamics_module).error_message_)

        output_value = pvecthermo[index_th]
        return output_value

    cpdef angular_distance(self, double z):
        """
        angular_distance(z)

        Return the angular diameter distance (exactly, the quantity defined by Class
        as index_bg_ang_distance in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_ang_distance
        background_module = deref(self._thisptr).GetBackgroundModule()
        index_bg_ang_distance = deref(background_module).index_bg_ang_distance_
        return self.get_background_value_at_z(z, index_bg_ang_distance)

    cpdef scale_independent_growth_factor(self, z):
        """
        scale_independent_growth_factor(z)

        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_D
        background_module = deref(self._thisptr).GetBackgroundModule()
        index_bg_D = deref(background_module).index_bg_D_
        return self.get_background_value_at_z(z, index_bg_D)

    cpdef scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_f
        background_module = deref(self._thisptr).GetBackgroundModule()
        index_bg_f = deref(background_module).index_bg_f_
        return self.get_background_value_at_z(z, index_bg_f)

    @cython.cdivision(True)
    cpdef z_of_tau(self, double tau):
        """
        Redshift corresponding to a given conformal time.

        Parameters
        ----------
        tau : float
                Conformal time
        """
        cdef:
            double output_value
            double z
            int bg_size_short
            int index_bg_a
            int last_index
            int status
            short inter_normal
            short short_info
            vector[double] pvecback

        short_info = self.ba.short_info
        inter_normal = self.ba.inter_normal

        background_module = deref(self._thisptr).GetBackgroundModule()
        bg_size_short = deref(background_module).bg_size_short_
        index_bg_a = deref(background_module).index_bg_a_

        pvecback.resize(bg_size_short)
        status = deref(background_module).background_at_tau(tau, short_info, inter_normal, &last_index, &pvecback[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)
        z = 1./pvecback[index_bg_a] - 1.
        return z

    cpdef Hubble(self, double z):
        """
        Hubble(z)

        Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_H
        background_module = deref(self._thisptr).GetBackgroundModule()
        index_bg_H = deref(background_module).index_bg_H_
        return self.get_background_value_at_z(z, index_bg_H)

    cpdef Om_m(self, double z):
        """
        Omega_m(z)

        Return the matter density fraction (exactly, the quantity defined by Class as index_bg_Omega_m
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_bg_Omega_m
        background_module = deref(self._thisptr).GetBackgroundModule()
        index_bg_Omega_m = deref(background_module).index_bg_Omega_m_
        return self.get_background_value_at_z(z, index_bg_Omega_m)


    cpdef ionization_fraction(self, double z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_th_xe
        thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
        index_th_xe = deref(thermodynamics_module).index_th_xe_
        return self.get_thermodynamics_value_at_z(z, index_th_xe)

    cpdef baryon_temperature(self, double z):
        """
        baryon_temperature(z)

        Give the baryon temperature for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef int index_th_Tb
        thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
        index_th_Tb = deref(thermodynamics_module).index_th_Tb_
        return self.get_thermodynamics_value_at_z(z, index_th_Tb)

    cpdef T_cmb(self):
        """
        Return the CMB temperature
        """
        return self.ba.T_cmb

    # redundent with a previous Omega_m() funciton,
    # but we leave it not to break compatibility
    cpdef Omega0_m(self):
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        return self.Omega_m()

    cpdef get_background(self):
        """
        Return an array of the background quantities at all times.

        Parameters
        ----------

        Returns
        -------
        background : dictionary containing background.
        """
        cdef:
            dict background
            double[:,::1] view
            int i
            int number_of_titles
            int status
            int timesteps
            string titles
            vector[double] data

        background_module = deref(self._thisptr).GetBackgroundModule()
        titles.resize(_MAXTITLESTRINGLENGTH_)
        status = deref(background_module).background_output_titles(<char*> titles.c_str())
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)

        tmp = <bytes> titles.c_str()
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(background_module).bt_size_

        data.resize(timesteps*number_of_titles)
        status = deref(background_module).background_output_data(number_of_titles, &data[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(background_module).error_message_)

        background = {}
        view = <double[:timesteps,:number_of_titles]> &data[0]
        for i in range(number_of_titles):
            background[names[i]] = np.asarray(view[:, i]).copy()

        return background

    cpdef get_thermodynamics(self):
        """
        Return the thermodynamics quantities.

        Returns
        -------
        thermodynamics : dictionary containing thermodynamics.
        """
        cdef:
            dict thermodynamics
            double[:,::1] view
            int i
            int number_of_titles
            int status
            int timesteps
            string titles
            vector[double] data

        thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()

        titles.resize(_MAXTITLESTRINGLENGTH_)
        status = deref(thermodynamics_module).thermodynamics_output_titles(<char*> titles.c_str())
        if status == _FAILURE_:
            raise CosmoSevereError(deref(thermodynamics_module).error_message_)

        tmp = <bytes> titles.c_str()
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(thermodynamics_module).tt_size_

        data.resize(timesteps*number_of_titles)

        status = deref(thermodynamics_module).thermodynamics_output_data(number_of_titles, &data[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(thermodynamics_module).error_message_)

        thermodynamics = {}
        view = <double[:timesteps,:number_of_titles]> &data[0]
        for i in range(number_of_titles):
            thermodynamics[names[i]] = np.asarray(view[:, i]).copy()

        return thermodynamics

    cpdef get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """
        cdef:
            string titles
            vector[double] data
            int status
            int timesteps
            int number_of_titles
            double[:,::1] view
            dict primordial
            int i

        primordial_module = deref(self._thisptr).GetPrimordialModule()

        titles.resize(_MAXTITLESTRINGLENGTH_)

        status = deref(primordial_module).primordial_output_titles(<char*> titles.c_str())
        if status == _FAILURE_:
            raise CosmoSevereError(deref(primordial_module).error_message_)

        tmp = <bytes> titles.c_str()
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(primordial_module).lnk_size_

        data.resize(timesteps*number_of_titles)
        status = deref(primordial_module).primordial_output_data(number_of_titles, &data[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(primordial_module).error_message_)

        primordial = {}
        view = <double[:timesteps,:number_of_titles]> &data[0]
        for i in range(number_of_titles):
            primordial[names[i]] = np.asarray(view[:, i]).copy()

        return primordial


    @cython.returns(dict)
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.ccall
    cpdef get_perturbations(self):
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
            Py_ssize_t i
            Py_ssize_t j
            Py_ssize_t k_size
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            dict tmpdict
            double** thedata
            double[:,::1] view
            int* thesizes
            list names
            list tmparray

        perturbation_module = deref(self._thisptr).GetPerturbationsModule()
        k_size = self.pt.k_output_values_num
        for mode in ['scalar','vector','tensor']:
            if mode=='scalar' and self.pt.has_scalars:
                thetitles = <bytes> deref(perturbation_module).scalar_titles_
                thedata = <double**> deref(perturbation_module).scalar_perturbations_data_
                thesizes = <int*> deref(perturbation_module).size_scalar_perturbation_data_
            elif mode=='vector' and self.pt.has_vectors:
                thetitles = <bytes> deref(perturbation_module).vector_titles_
                thedata = <double**> deref(perturbation_module).vector_perturbations_data_
                thesizes = <int*> deref(perturbation_module).size_vector_perturbation_data_
            elif mode=='tensor' and self.pt.has_tensors:
                thetitles = <bytes> deref(perturbation_module).tensor_titles_
                thedata = <double**> deref(perturbation_module).tensor_perturbations_data_
                thesizes = <int*> deref(perturbation_module).size_tensor_perturbation_data_
            else:
                continue
            thetitles = str(thetitles.decode())
            names = thetitles.strip("\t").split("\t")
            number_of_titles = len(names)
            tmparray = []
            if number_of_titles == 0:
                continue
            for j in range(k_size):
                timesteps = thesizes[j]//number_of_titles
                tmpdict = {}
                view = <double[:timesteps,:number_of_titles]> thedata[j]
                for i in range(number_of_titles):
                    tmpdict[names[i]] = np.asarray(view[:,i])
                tmparray.append(tmpdict)
            perturbations[mode] = tmparray

        return perturbations

    cpdef get_transfer(self, double z=0., output_format='class'):
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
        cdef:
            string titles
            vector[double] data
            char ic_info[1024]
            FileName ic_suffix
            file_format outf
            int status
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            Py_ssize_t size_ic_data
            Py_ssize_t ic_num
            dict tmpdict
            dict transfers
            double[:,:,::1] view
            int index_ic
            int i

        if (not self.pt.has_density_transfers) and (not self.pt.has_velocity_transfers):
            return {}

        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        perturbations_module = deref(self._thisptr).GetPerturbationsModule()
        index_md = deref(perturbations_module).index_md_scalars_
        titles.resize(_MAXTITLESTRINGLENGTH_)

        status = deref(perturbations_module).perturb_output_titles(outf, <char*> titles.c_str())
        if status == _FAILURE_:
            raise CosmoSevereError(deref(perturbations_module).error_message_)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.strip("\t").split("\t")
        number_of_titles = len(names)
        timesteps = deref(perturbations_module).k_size_[index_md]

        size_ic_data = timesteps*number_of_titles
        ic_num = deref(perturbations_module).ic_size_[index_md]

        data.resize(size_ic_data*ic_num)

        status = deref(perturbations_module).perturb_output_data(outf, z, number_of_titles, &data[0])
        if status == _FAILURE_:
            raise CosmoSevereError(deref(perturbations_module).error_message_)

        transfers = {}
        view = <double[:ic_num,:timesteps,:number_of_titles]> &data[0]

        for index_ic in range(ic_num):
            status = deref(perturbations_module).perturb_output_firstline_and_ic_suffix(index_ic, ic_info, ic_suffix)
            if status == _FAILURE_:
                raise CosmoSevereError(deref(perturbations_module).error_message_)
            ic_key = <bytes> ic_suffix
            ic_key = str(ic_key.decode())

            tmpdict = {}
            for i in range(number_of_titles):
                tmpdict[names[i]] = np.asarray(view[index_ic, :, i]).copy()

            if ic_num == 1:
                transfers = tmpdict
            else:
                transfers[ic_key] = tmpdict

        return transfers

    @cython.cdivision(True)
    cdef get_slowroll_parameters(self):
        cdef:
            (double, double, double) epsilons
            double C
            double alpha_s
            double eps1
            double eps2
            double eps23
            double n_s
            double r

        primordial_module = deref(self._thisptr).GetPrimordialModule()
        r = deref(primordial_module).r_
        n_s = deref(primordial_module).n_s_
        alpha_s = deref(primordial_module).alpha_s_
        C = EULER + LOGE2 - 2.0

        eps1 = r*(1./16. + C/16.*(r/8. + n_s - 1.))
        eps2 = -n_s + 1. + C*alpha_s - r*(1./8. + 1./8.*(n_s - 1.)*(C - 1.5)) - (r/8.)**2*(C - 1.)
        eps23 = 1./8.*(r**2/8. + (n_s - 1.)*r - 8.*alpha_s)
        epsilons = eps1, eps2, eps23
        return epsilons

    @cython.cdivision(True)
    cpdef get_current_derived_parameters(self, names):
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
        cdef:
            dict derived
            double value
            double eps1
            double eps2
            double eps23
            (double, double, double) epsilons

        if not isinstance(names, list):
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
                value = self.age()
            elif name == 'conformal_age':
                background_module = deref(self._thisptr).GetBackgroundModule()
                value = deref(background_module).conformal_age_
            #TODO: Need to wrap NCDM class to get the masses
            #elif name == 'm_ncdm_in_eV':
            #    value = self.ba.m_ncdm_in_eV[0]
            elif name == 'm_ncdm_tot':
                # TODO: Is this really what we want??
                value = self.ba.Omega0_ncdm_tot*self.ba.h*self.ba.h*93.14
            elif name == 'Neff':
                value = self.Neff()
            elif name == 'Omega_m':
                value = self.Omega_m()
            elif name == 'omega_m':
                value = self.Omega_m()/self.ba.h/self.ba.h
            elif name == 'xi_idr':
                value = self.ba.T_idr/self.ba.T_cmb
            elif name == 'N_dg':
                value = self.ba.Omega0_idr/self.ba.Omega0_g*8./7.*pow(11./4.,4./3.)
            elif name == 'Gamma_0_nadm':
                value = self.th.a_idm_dr*(4./3.)*(self.ba.h*self.ba.h*self.ba.Omega0_idr)
            elif name == 'a_dark':
                value = self.th.a_idm_dr
            elif name == 'tau_reio':
                value = self.tau_reio()
            elif name == 'z_reio':
                value = self.z_reio()
            elif name == 'z_rec':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).z_rec_
            elif name == 'tau_rec':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).tau_rec_
            elif name == 'rs_rec':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).rs_rec_
            elif name == 'rs_rec_h':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).rs_rec_*self.ba.h
            elif name == 'ds_rec':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ds_rec_
            elif name == 'ds_rec_h':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ds_rec_*self.ba.h
            elif name == 'ra_rec':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ra_rec_
            elif name == 'ra_rec_h':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ra_rec_*self.ba.h
            elif name == 'da_rec':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).da_rec_
            elif name == 'da_rec_h':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).da_rec_*self.ba.h
            elif name == 'z_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).z_star_
            elif name == 'tau_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).tau_star_
            elif name == 'rs_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).rs_star_
            elif name == 'ds_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ds_star_
            elif name == 'ra_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ra_star_
            elif name == 'da_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).da_star_
            elif name == 'rd_star':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).rd_star_
            elif name == 'z_d':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).z_d_
            elif name == 'tau_d':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).tau_d_
            elif name == 'ds_d':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ds_d_
            elif name == 'ds_d_h':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).ds_d_*self.ba.h
            elif name == 'rs_d':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).rs_d_
            elif name == 'rs_d_h':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                rs_d = deref(thermodynamics_module).rs_d_
                value = rs_d*self.ba.h
            elif name == '100*theta_s':
                value = self.theta_s_100()
            elif name == '100*theta_star':
                value = self.theta_star_100()
            elif name == 'YHe':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).YHe_
            elif name == 'n_e':
                thermodynamics_module = deref(self._thisptr).GetThermodynamicsModule()
                value = deref(thermodynamics_module).n_e_
            elif name == 'A_s':
                value = self.A_s()
            elif name == 'ln10^{10}A_s':
                value = log(1.e10*self.A_s())
            elif name == 'n_s':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).n_s_
            elif name == 'alpha_s':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).alpha_s_
            elif name == 'beta_s':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).beta_s_
            elif name == 'r':
                # This is at the pivot scale
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).r_
            elif name == 'r_0002':
                # at k_pivot = 0.002/Mpc
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                r = deref(primordial_module).r_
                n_s = deref(primordial_module).n_s_
                n_t = deref(primordial_module).n_t_
                alpha_s = deref(primordial_module).alpha_s_
                value = r*(0.002/self.pm.k_pivot)**(n_t - n_s - 1 + 0.5*alpha_s*log(0.002/self.pm.k_pivot))
            elif name == 'n_t':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).n_t_
            elif name == 'alpha_t':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).alpha_t_
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
                epsilons = self.get_slowroll_parameters()
                eps1, eps2, eps23 = epsilons
                value = eps1*((1. - eps1/3. + eps2/6.)/(1. - eps1/3.))**2
            elif name == 'eta_V':
                epsilons = self.get_slowroll_parameters()
                eps1, eps2, eps23 = epsilons
                value = (2.*eps1 - eps2/2. - 2./3.*eps1**2 + 5./6.*eps1*eps2 - eps2**2/12. - eps23/6.)/(1. - eps1/3.)
            elif name == 'ksi_V^2':
                epsilons = self.get_slowroll_parameters()
                eps1, eps2, eps23 = epsilons
                value = 2.*(1. - eps1/3. + eps2/6.)*(2.*eps1**2 - 3./2.*eps1*eps2 + eps23/4.)/(1. - eps1/3.)**2
            elif name == 'exp_m_2_tau_As':
                value = exp(-2*self.tau_reio())*self.A_s()
            elif name == 'phi_min':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).phi_min_
            elif name == 'phi_max':
                primordial_module = deref(self._thisptr).GetPrimordialModule()
                value = deref(primordial_module).phi_max_
            elif name == 'sigma8':
                value = self.sigma8()
            elif name == 'sigma8_cb':
                value = self.sigma8_cb()
            elif name == 'k_eq':
                value = self.k_eq()
            else:
                raise CosmoSevereError("%s was not recognized as a derived parameter" % name)
            derived[name] = value
        return derived

    cpdef nonlinear_scale(self, double[::1] z, Py_ssize_t z_size):
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
        cdef:
            Py_ssize_t index_z
            int status
            double k_nl_val
            double k_nl_cb_val
            double[:] k_nl

        k_nl_arr = np.empty(z_size, np.double)
        k_nl = k_nl_arr
        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        for index_z in range(z_size):
            status = deref(nonlinear_module).nonlinear_k_nl_at_z(z[index_z], &k_nl_val, &k_nl_cb_val)
            if status == _FAILURE_:
                raise CosmoSevereError(deref(nonlinear_module).error_message_)
            k_nl[index_z] = k_nl_val

        return k_nl_arr

    cpdef nonlinear_scale_cb(self, double[::1] z, Py_ssize_t z_size):
        """

        nonlinear_scale_cb(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size

        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef:
            Py_ssize_t index_z
            int status
            double k_nl_val
            double k_nl_cb_val
            double[:] k_nl_cb

        k_nl_cb_arr = np.empty(z_size, np.double)
        k_nl_cb = k_nl_cb_arr
        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        for index_z in range(z_size):
            status = deref(nonlinear_module).nonlinear_k_nl_at_z(z[index_z], &k_nl_val, &k_nl_cb_val)
            if status == _FAILURE_:
                raise CosmoSevereError(deref(nonlinear_module).error_message_)
            k_nl_cb[index_z] = k_nl_cb_val

        return k_nl_cb_arr

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

        # Set the module to the current values
        self._pars = data.cosmo_arguments
        self.reset()
        self.compute()

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

    cdef get_pk_array_general(self, double[::1] k,  double[::1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef:
            cdef double[::1] pk
            cdef double[::1] pk_cb
            int status
            pk_outputs linear_or_nonlinear
        if nonlinear == 0:
            linear_or_nonlinear = pk_linear
        else:
            linear_or_nonlinear = pk_nonlinear

        pk_arr = np.empty(k_size*z_size, np.double)
        pk = pk_arr
        pk_cb_arr = np.empty(k_size*z_size, np.double)
        pk_cb = pk_cb_arr

        nonlinear_module = deref(self._thisptr).GetNonlinearModule()
        deref(nonlinear_module).nonlinear_pks_at_kvec_and_zvec(pk_linear, &k[0], k_size, &z[0], z_size, &pk[0],  &pk_cb[0])

        return pk_arr, pk_cb_arr

    cpdef get_pk_array(self, double[::1] k,  double[::1] z, int k_size, int z_size, nonlinear):
        return self.get_pk_array_general(k, z, k_size, z_size, nonlinear)[0]

    cpdef get_pk_cb_array(self, double[::1] k,  double[::1] z, int k_size, int z_size, nonlinear):
        return self.get_pk_array_general(k, z, k_size, z_size, nonlinear)[1]

    cpdef Omega0_k(self):
        """ Curvature contribution """
        return self.ba.Omega0_k

    cpdef Omega0_cdm(self):
        return self.ba.Omega0_cdm
