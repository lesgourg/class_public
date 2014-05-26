"""
.. module:: classy
    :synopsis: Python wrapper around CLASS
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>

This module defines a class called Class. It is used with Monte Python to
extract cosmological parameters.

"""
from math import exp,log
import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
cimport cython

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i

# Import the .pxd containing definitions
from cclassy cimport *

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=""):
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
    pass


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
    cdef file_content fc

    cpdef int ready # Flag
    cpdef object _pars # Dictionary of the parameters
    cpdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

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
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()

        # TEST
        #raise CosmoSevereError

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
        self.ready=False

    def cleanup(self):
        if self.ready==False:
            return True
        for i in range(len(self._pars)):
            if self.fc.read[i]==0:
                del(self._pars[self.fc.name[i]])

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

            dumc = kk
            sprintf(self.fc.name[i],"%s",dumc)
            dumcp = str(self._pars[kk])
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = _FALSE_
            i+=1

    # Called at the end of a run, to free memory
    def struct_cleanup(self):
        if self.ready == _FALSE_:
             return
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

    # Ensure the full module dependency
    def _check_task_dependency(self,lvl):
        if "lensing" in lvl:
            lvl.append("spectra")
        if "spectra" in lvl:
            lvl.append("transfer")
        if "transfer" in lvl:
            lvl.append("nonlinear")
        if "nonlinear" in lvl:
            lvl.append("primordial")
        if "primordial" in lvl:
            lvl.append("perturb")
        if "perturb" in lvl:
            lvl.append("thermodynamics")
        if "thermodynamics" in lvl:
            lvl.append("background")
        if len(lvl)!=0 :
            lvl.append("input")
        return lvl

    def _pars_check(self,key,value,contains=False,add=""):
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

    # Main function, computes (call _init methods for all desired modules). This
    # is called in MontePython, and this ensures that the Class instance of this
    # class contains all the relevant quantities. Then, one can deduce Pk, Cl,
    # etc...
    # lvl default value should be left as an array (was creating problem when
    # casting as a set later on in check_task_dependency) (TOFIX more)
    def compute(self, lvl=["lensing"]):
        """
        compute(lvl=["lensing"])

        Compute CMB quantities up to module lvl

        Parameters
        ----------
        lvl : list
                last module desired. For instance, "lensing". Then Class will compute all
                the modules required to compute this desired module through the
                (internal) function _check_task_dependency

        """
        cdef ErrorMsg errmsg
        cdef int ierr
        cdef char* dumc

        lvl = self._check_task_dependency(lvl)

        if self.ready and self.ncp.issuperset(lvl):
            return

        self.ready = False

        self._fillparfile()

        # empty all
        self.ncp=set()

        # compute
        if "input" in lvl:
            ierr = input_init(
                &self.fc,
                &self.pr,
                &self.ba,
                &self.th,
                &self.pt,
                &self.tr,
                &self.pm,
                &self.sp,
                &self.nl,
                &self.le,
                &self.op,
                errmsg)
            if ierr == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.ncp.add("input")

            problem_flag = False
            problematic_parameters = []
            for i in range(self.fc.size):
                if self.fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(self.fc.name[i])
            if problem_flag:
                raise CosmoSevereError(
                    "Class did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))

        if "background" in lvl:
            if background_init(&(self.pr),&(self.ba)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.ba.error_message)
            self.ncp.add("background")

        if "thermodynamics" in lvl:
            if thermodynamics_init(&(self.pr),&(self.ba),&(self.th)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.th.error_message)
            self.ncp.add("thermodynamics")

        if "perturb" in lvl:
            if perturb_init(&(self.pr),&(self.ba),&(self.th),&(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            self.ncp.add("perturb")

        if "primordial" in lvl:
            if primordial_init(&(self.pr),&(self.pt),&(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            self.ncp.add("primordial")

        if "nonlinear" in lvl:
            if nonlinear_init(&self.pr,&self.ba,&self.th,&self.pt,&self.pm,&self.nl) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.nl.error_message)
            self.ncp.add("nonlinear")

        if "transfer" in lvl:
            if transfer_init(&(self.pr),&(self.ba),&(self.th),&(self.pt),&(self.nl),&(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)
            self.ncp.add("transfer")

        if "spectra" in lvl:
            if spectra_init(&(self.pr),&(self.ba),&(self.pt),&(self.pm),&(self.nl),&(self.tr),&(self.sp)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sp.error_message)
            self.ncp.add("spectra")

        if "lensing" in lvl:
            if lensing_init(&(self.pr),&(self.pt),&(self.sp),&(self.nl),&(self.le)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.le.error_message)
            self.ncp.add("lensing")

        self.ready = True

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
        cdef double *rcl = <double*> calloc(self.sp.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md[index_md] = <double*> calloc(self.sp.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.sp.ct_size*self.sp.ic_ic_size[index_md], sizeof(double))

        lmaxR = self.sp.l_max_tot
        if lmax==-1:
            lmax=lmaxR
        if lmax>lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}
        for elem in ['tt','te','ee','bb','pp','tp']:
            cl[elem] = np.ndarray(lmax+1, dtype=np.double)
            cl[elem][:2]=0
        for ell from 2<=ell<lmax+1:
            if spectra_cl_at_l(&self.sp,ell,rcl,cl_md,cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.sp.error_message)
            cl['tt'][ell] = rcl[self.sp.index_ct_tt]
            cl['te'][ell] = rcl[self.sp.index_ct_te]
            cl['ee'][ell] = rcl[self.sp.index_ct_ee]
            cl['bb'][ell] = rcl[self.sp.index_ct_bb]
            cl['pp'][ell] = rcl[self.sp.index_ct_pp]
            cl['tp'][ell] = rcl[self.sp.index_ct_tp]

        free(rcl)
        return cl

    def lensed_cl(self, lmax=-1,nofail=False):
        """
        lensed_cl(lmax=-1, nofail=False)

        Return a dictionary of the lensed C_l

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
        lmaxR = self.le.l_lensed_max

        if lmax==-1:
            lmax=lmaxR
        if lmax>lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}
        for elem in ['tt','te','ee','bb','pp','tp']:
            cl[elem] = np.ndarray(lmax+1, dtype=np.double)
            cl[elem][:2]=0
        for ell from 2<=ell<lmax+1:
            if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
                raise CosmoSevereError(self.le.error_message)
            cl['tt'][ell] = lcl[self.le.index_lt_tt]
            cl['te'][ell] = lcl[self.le.index_lt_te]
            cl['ee'][ell] = lcl[self.le.index_lt_ee]
            cl['bb'][ell] = lcl[self.le.index_lt_bb]
            cl['pp'][ell] = lcl[self.le.index_lt_pp]
            cl['tp'][ell] = lcl[self.le.index_lt_tp]

        free(lcl)
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

    # Gives the pk for a given (k,z)
    def pk(self,double k,double z):
        """
        Gives the pk for a given k and z

        .. note::

            there is an additional check to verify if output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk
        cdef double pk_velo
        cdef double pk_cross
        cdef int dummy

        # Quantities for the isocurvature modes
        cdef double *pk_ic = <double*> calloc(self.sp.ic_ic_size[self.sp.index_md_scalars], sizeof(double))
        abort = True
        if 'output' in self._pars:
            options = self._pars['output'].split()
            for option in options:
                if option in ['mPk', 'mTk', 'vTk']:
                    abort = False
                    break
        if abort:
            raise CosmoSevereError(
                "No power spectrum nor transfer function"
                " asked: you should not ask for a power"
                " spectrum, then")

        if (self.nl.method == 0):
             if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,pk_ic)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)
        else:
             if spectra_pk_nl_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk) ==_FAILURE_:
                    raise CosmoSevereError(self.sp.error_message)
        #free(junk)
        return pk

    def get_pk(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        #cdef np.ndarray[DTYPE_t, ndim = 3] k_int = k
        #cdef np.ndarray[DTYPE_t, ndim = 1] z_int = z
        #cdef np.ndarray k_int
        #cdef np.ndarray z_int

        #k_int = np.zeros((k_size,z_size,mu_size),'float64')
        #z_int = np.zeros((z_size),'float64')
        #k_int = k
        #z_int = z
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self._pk(k[index_k,index_z,index_mu],z[index_z])
        return pk

    # Avoids using hardcoded numbers for tt, te, ... indexes in the tables.
    def return_index(self):
        index = {}
        index['tt'] = self.le.index_lt_tt
        index['te'] = self.le.index_lt_te
        index['ee'] = self.le.index_lt_ee
        index['bb'] = self.le.index_lt_bb
        index['pp'] = self.le.index_lt_pp
        index['tp'] = self.le.index_lt_tp
        return index

    def age(self):
        self.compute(["background"])
        return self.ba.age

    def h(self):
        return self.ba.h

    def n_s(self):
        return self.pm.n_s

    # Defined twice ?
    def Omega_m(self):
        return self.ba.Omega0_b+self.ba.Omega0_cdm

    def Omega_b(self):
        return self.ba.Omega0_b

    def omega_b(self):
        return self.ba.Omega0_b * self.ba.h * self.ba.h

    def Neff(self):
        return self.ba.Neff

    def sigma8(self):
        self.compute(["spectra"])
        return self.sp.sigma8

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
        Return the sum of Omega0 for baryon and CDM
        """
        return self.ba.Omega0_b+self.ba.Omega0_cdm

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
        #cdef double *k_nl

        #k_nl = <double*> calloc(z_size,sizeof(double))
        for index_z in range(z_size):
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return k_nl

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
        self.get_current_derived_parameters(data)
        for elem in data.get_mcmc_parameters(['derived']):
            data.mcmc_parameters[elem]['current'] /= \
                data.mcmc_parameters[elem]['scale']
            params[elem] = data.mcmc_parameters[elem]['current']

        ctx.add('boundary', True)
        # Store itself into the context, to be accessed by the likelihoods
        ctx.add('cosmo', self)
