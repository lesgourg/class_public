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
from math import exp,log
import numpy as np
from os.path import abspath, dirname
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

# Nils : Added for python 3.x and python 2.x compatibility
import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

ctypedef np.float64_t DTYPE_t
ctypedef np.int32_t DTYPE_i

# A little utility function to ensure backwards compatibility. This allows classy to define properties that you can also call at the same time
class CallableFloat(float):
  def __call__(self):
    return self


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

    cdef int computed # Flag to see if classy has already computed with the given pars
    cdef int allocated # Flag to see if classy structs are allocated already
    cdef object _pars # Dictionary of the parameters
    cdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

    cdef char path_to_this[1000]

    _levellist = ["input","background","thermodynamics","perturbations", "primordial", "fourier", "transfer", "harmonic", "lensing", "distortions"]

    #############################################################################
    # Below we define the list of all properties to be used in CLASS            #
    # They can either be used as a property (e.g. cosmo.h)                      #
    # or as a function (cosmo.h() ) in order to keep backwards compatibility    #
    #############################################################################

    # Energy density fractional properties
    @property
    def Omega_nu(self):
      """Return the fractional massive neutrino/warm dark matter (non-cold dark matter) density today"""
      return CallableFloat(self.ba.Omega0_ncdm_tot)
    @property
    def Omega_m(self):
      """Return the fractional non-relativistic matter density today"""
      return CallableFloat(self.ba.Omega0_m)
    @property
    def Omega_r(self):
      """Return the fractional radiation density today"""
      return CallableFloat(self.ba.Omega0_r)
    @property
    def Omega_Lambda(self):
      """Return the effective fractional contribution of the cosmological constant to the density today"""
      return CallableFloat(self.ba.Omega0_lambda)
    @property
    def Omega_g(self):
      """Return the fractional photon density today"""
      return CallableFloat(self.ba.Omega0_g)
    @property
    def Omega_b(self):
      """Return the fractional baryon density today"""
      return CallableFloat(self.ba.Omega0_b)
    @property
    def omega_b(self):
      """Return the reduced fractional baryon density today"""
      return CallableFloat(self.ba.Omega0_b*self.ba.h*self.ba.h)
    @property
    def Omega_k(self):
      """Return the effective fractional contribution of curvature to the density today"""
      return CallableFloat(self.ba.Omega0_k)
    @property
    def Omega_cdm(self):
      """Return the fractional cold dark matter density today"""
      return CallableFloat(self.ba.Omega0_cdm)
    @property
    def Omega_fld(self):
      """Return the fractional fluid density today"""
      return CallableFloat(self.ba.Omega0_fld)

    # Other properties related to the background
    @property
    def H0(self):
      """Return the Hubble constant today in units of km/s/Mpc"""
      return CallableFloat(self.ba.H0)
    @property
    def h(self):
      """Return the reduced Hubble constant today (=H0 in units of 100km/s/Mpc)"""
      return CallableFloat(self.ba.h)
    @property
    def age(self):
      """Return the age of the universe in Gigayears"""
      self.compute(["background"])
      return CallableFloat(self.ba.age)
    @property
    def z_eq(self):
      """Return the redshift of the time of equality between radiation and matter (dimensionless)"""
      self.compute(["background"])
      return CallableFloat(1./self.ba.a_eq-1.)
    @property
    def k_eq(self):
      """Return the wavenumber such that k=aH at the time of equality between radiation and matter (units of 1/Mpc)"""
      self.compute(["background"])
      return CallableFloat(self.ba.a_eq*self.ba.H_eq)
    @property
    def Neff(self):
      """Return the effective neutrino number N_eff (dimensionless) which parametrizes the density of radiation in the early universe, before non-cold dark matter particles become non-relativistic. Should be 3.044 in the standard model."""
      self.compute(["background"])
      return CallableFloat(self.ba.Neff)
    @property
    def T_cmb(self):
      """Return the photon temperature T_cmb (units of K) evaluated at z=0"""
      return CallableFloat(self.ba.T_cmb)

    # Other properties related to thermodynamical quantities
    @property
    def tau_reio(self):
      """Return the optical depth of reionization, which is the
         average number of re-scatterings any single CMB photon
         is subject to during re-ionization."""
      self.compute(["thermodynamics"])
      return CallableFloat(self.th.tau_reio)
    @property
    def z_reio(self):
      """Return the reionization redshift (dimensionless). This
         redshift is a free parameter in the analytic form assummed for
         the free electron fraction x_e(z). It is either passed by the
         user, or adjusted by CLASS using a shooting method in order to
         match a required value of the reionization optical depth."""
      self.compute(["thermodynamics"])
      return CallableFloat(self.th.z_reio)
    @property
    def rs_drag(self):
      """Return the comoving sound horizon (units of Mpc) at baryon
        drag time. CLASS defines this time as the time when the baryon
        optical depth crosses 1. This baryon optical depth is obtained
        by integrating the Thomson scattering rate of baryons, that
        is, the Thomson scattering rate of photons rescaled by 1/R,
        where R = 3 rho_b / 4 rho_gamma."""
      self.compute(["thermodynamics"])
      return CallableFloat(self.th.rs_d)
    @property
    def theta_s_100(self):
      """Return the angular scale of the sound horizon at recombination
        multiplied by 100 (in radians, that is, dimensionless). The
        function uses theta_s = d_s(z_rec)/d_a(t_rec) =
        r_s(z_rec)/r_a(z_rec), where the d are physical distances and
        the r are comoving distances. CLASS defines recombination as
        the time at which the photon visibility function reaches its maximum."""
      self.compute(["thermodynamics"])
      return CallableFloat(100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec))
    @property
    def theta_star_100(self):
      """Return the angular scale of the sound horizon at photon
        decoupling multiplied by 100 (in radian, that is,
        dimensionless). The function uses theta_s =
        d_s(z_star)/d_a(t_star) = r_s(z_star)/r_a(z_star), where the d
        are physical distances and the r are comoving distances. CLASS
        defines decoupling as the time at which the photon optical
        depth crosses one."""
      self.compute(["thermodynamics"])
      return CallableFloat(100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star))

    # Other properties related to the peruturbations
    @property
    def n_s(self):
      """Return the scalar tilt n_s (dimensionless) of the primordial scalar spectrum"""
      return CallableFloat(self.pm.n_s)
    @property
    def simga8(self):
      """Return sigma8 (dimensionless), the root mean square (rms) of the relative density fluctuation of
        total matter in spheres of radius R= 8 h/Mpc at z=0.
        To compute this, one needs that the 'ouput' field contains at least 'mPk'."""
      self.compute(["fourier"])
      if (self.pt.has_pk_matter == _FALSE_):
          raise CosmoSevereError("No power spectrum computed. In order to get sigma8, you must add mPk to the list of outputs.")
      return CallableFloat(self.fo.sigma8[self.fo.index_pk_m])
    @property
    def simga8_cb(self):
      """Return sigma8_cb (dimensionless), the root mean square (rms) of the relative density fluctuation of
        CDM+baryons in spheres of radius R= 8 h/Mpc at z=0
        To compute this, one needs that the 'ouput' field contains at least 'mPk'."""
      self.compute(["fourier"])
      if (self.pt.has_pk_matter == _FALSE_):
          raise CosmoSevereError("No power spectrum computed. In order to get sigma8_cb, you must add mPk to the list of outputs.")
      return CallableFloat(self.fo.sigma8[self.fo.index_pk_cb])
    @property
    def S8(self):
      """Return S8 = sigma8*(Omega_m/0.3)^0.5 (dimensionless), the root mean square (rms) of the relative density fluctuation of
        total matter in spheres of radius R= 8 h/Mpc at z=0 multiplied by (Omega_m/0.3)^0.5"""
      self.compute(["fourier"])
      return CallableFloat(self.sigma8()*np.sqrt(self.Omega_m()/0.3))
    @property
    def nonlinear_method(self):
      """Return which method is being used to compute non-linear corrections (as an integer, see the corresponding enum in the C code to check for the correspondence."""
      return CallableFloat(self.fo.method)

    # Unit convrsion proprtis
    @property
    def density_factor(self):
        """Returns the factor to convert from the CLASS units of density, (8piG/3c^2 rho) in Mpc^-2,
        to the actual density rho in kg/m^3 (SI units)"""
        return 3*_c_*_c_/(8*np.pi*_G_)/(_Mpc_over_m_*_Mpc_over_m_)
    @property
    def Mpc_to_m(self):
        """Returns the factor to convert from Megaparsecs (Mpc) to meters (m)
        (factor defined internally in CLASS)"""
        return _Mpc_over_m_
    @property
    def kg_to_eV(self):
        """Returns the factor to convert from kilograms (kg) to electron-volts (eV) when using natural units (inferred from factors defined internally in CLASS)"""
        return _c_*_c_/_eV_
    @property
    def c(self):
        """Returns the speed of light in m/s (factor defined internally in CLASS)"""
        return _c_
    @property
    def kgm3_to_eVMpc3(self):
        """Returns the factor to convert from kilograms per meter cube (kg/m^3) to electron-volts per Megaparsec cube (eV/Mpc^3) when using natural units"""
        return self.kg_to_eV*self.Mpc_to_m**3
    @property
    def kg_to_Msol(self):
        """Returns the factor to convert from kilograms to solar masses"""
        return 1/(2.0e30)
    @property
    def kgm3_to_MsolMpc3(self):
        """Returns the factor to convert from kilograms per meter cube (kg/m^3) to solar mass per Megaparsec cube (Msol/Mpc^3)"""
        return self.kg_to_Msol*self.Mpc_to_m**3



    # Special properties
    @property
    def pars(self):
      return self._pars
    @property
    def state(self):
      return True

    ###############################################################
    # Now we can start with the actual code describing the classy #
    ###############################################################

    def set_default(self):
        _pars = {
            "output":"tCl mPk",}
        self.set(**_pars)

    def __cinit__(self, default=False):
        cdef char* dumc
        self.allocated = False
        self.computed = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()
        try:
          import importlib.resources
          resource_path = abspath(importlib.resources.files('classy'))
        except ImportError as ie:
          resource_path = dirname(abspath(__file__))
        path_to_this_as_bytes = resource_path.encode()
        dumc = path_to_this_as_bytes
        sprintf(self.path_to_this,"%s",dumc)

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

        self.fc.size = len(self._pars)+1 if 'base_path' not in self._pars else len(self._pars)
        self.fc.name = <FileArg*> malloc(sizeof(FileArg)*self.fc.size)
        assert(self.fc.name!=NULL)

        self.fc.value = <FileArg*> malloc(sizeof(FileArg)*self.fc.size)
        assert(self.fc.value!=NULL)

        self.fc.read = <short*> malloc(sizeof(short)*self.fc.size)
        assert(self.fc.read!=NULL)

        # fill parameter file
        i = 0
        for kk in self._pars:

            dumcp = kk.strip().encode()
            dumc = dumcp
            sprintf(self.fc.name[i],"%s",dumc)
            dumcp = str(self._pars[kk]).strip().encode()
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = _FALSE_
            i+=1
        if not 'base_path' in self._pars:
            dumcp = str('base_path').encode()
            dumc = dumcp
            sprintf(self.fc.name[i],"%s", dumc)
            dumc = self.path_to_this
            sprintf(self.fc.value[i], "%s", dumc)

    # Called at the end of a run, to free memory
    def struct_cleanup(self):
        if(self.allocated != True):
          return
        if self.sd.is_allocated:
            distortions_free(&self.sd)
        if self.le.is_allocated:
            lensing_free(&self.le)
        if self.hr.is_allocated:
            harmonic_free(&self.hr)
        if self.tr.is_allocated:
            transfer_free(&self.tr)
        if self.fo.is_allocated:
            fourier_free(&self.fo)
        if self.pm.is_allocated:
            primordial_free(&self.pm)
        if self.pt.is_allocated:
            perturbations_free(&self.pt)
        if self.th.is_allocated:
            thermodynamics_free(&self.th)
        if self.ba.is_allocated:
            background_free(&self.ba)
        self.ncp = set()

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
        # If it's a string only, treat as a list
        if isinstance(level, str):
          level=[level]
        # For each item in the list
        levelset = set()
        for item in level:
          # If the item is not in the list of allowed levels, make error message
          if item not in self._levellist:
            raise CosmoSevereError("Unknown computation level: '{}'".format(item))
          # Otherwise, add to list of levels up to and including the specified level
          levelset.update(self._levellist[:self._levellist.index(item)+1])
        return levelset

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

    def compute(self, level=["distortions"]):
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
            if input_read_from_file(&self.fc, &self.pr, &self.ba, &self.th,
                                    &self.pt, &self.tr, &self.pm, &self.hr,
                                    &self.fo, &self.le, &self.sd, &self.op, errmsg) == _FAILURE_:
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

        # The following list of computation is straightforward. If the "_init"
        # methods fail, call `struct_cleanup` and raise a CosmoComputationError
        # with the error message from the faulty module of CLASS.
        if "background" in level:
            if background_init(&(self.pr), &(self.ba)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.ba.error_message)
            self.ncp.add("background")

        if "thermodynamics" in level:
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.th.error_message)
            self.ncp.add("thermodynamics")

        if "perturbations" in level:
            if perturbations_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            self.ncp.add("perturbations")

        if "primordial" in level:
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            self.ncp.add("primordial")

        if "fourier" in level:
            if fourier_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.fo) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.fo.error_message)
            self.ncp.add("fourier")

        if "transfer" in level:
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.fo), &(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)
            self.ncp.add("transfer")

        if "harmonic" in level:
            if harmonic_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.fo), &(self.tr),
                            &(self.hr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.hr.error_message)
            self.ncp.add("harmonic")

        if "lensing" in level:
            if lensing_init(&(self.pr), &(self.pt), &(self.hr),
                            &(self.fo), &(self.le)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.le.error_message)
            self.ncp.add("lensing")

        if "distortions" in level:
            if distortions_init(&(self.pr), &(self.ba), &(self.th),
                                &(self.pt), &(self.pm), &(self.sd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sd.error_message)
            self.ncp.add("distortions")

        self.computed = True

        # At this point, the cosmological instance contains everything needed. The
        # following functions are only to output the desired numbers
        return

    def set_baseline(self, baseline_name):
        """
        Set all input parameters to predefined baseline settings

        This function offers an opportunity to quickly set all input parameters in one single line.

        Currently, the possible predefined setting names are:

        1. 'p18lb', or alternatively, any name containing 'planck', '18', 'lens', and 'bao'
           Sets all parameters for computing the lensed and unlensed CMB l's and the matter power spectrum P(k) at z=0
           (linear and non-linear with Halofit corrections)
           for the best-fit model of the case 'Planck TTTEEE+lens + BAO' of the Planck 2018 papers.

        2. 'p18l', or alternatively, any name containing 'planck', '18', and 'lens'
           Sets all parameters for computing the lensed and unlensed CMB l's and the matter power spectrum P(k) at z=0
           (linear and non-linear with Halofit corrections)
           for the best-fit model of the case 'Planck TTTEEE+lens' of the Planck 2018 papers.

        3. 'p18', or alternatively, any name containing 'planck' and '18'
           Sets all parameters for computing the lensed and unlensed CMB l's and the matter power spectrum P(k) at z=0
           (linear and non-linear with Halofit corrections)
           for the best-fit model of the case 'Planck TTTEEE' of the Planck 2018 papers.

        This choice of parameter settings is the same as in the baseline input files (input/baseline*.param)
        of montepython [https://github.com/brinckmann/montepython_public] (see also 1210.7183, 1804.07261)

        Parameters
        ----------
        baseline_name : str
            Predefined setting name

        Returns
        -------
        None
        """

        if ('planck' in baseline_name and '18' in baseline_name and 'lens' in baseline_name and 'bao' in baseline_name) or 'p18lb' in baseline_name.lower():
          self.set({'omega_b':2.255065e-02,
                    'omega_cdm':1.193524e-01,
                    'H0':6.776953e+01,
                    'A_s':2.123257e-09,
                    'n_s':9.686025e-01,
                    'z_reio':8.227371e+00,

                    'N_ur':2.0328,
                    'N_ncdm':1,
                    'm_ncdm':0.06,
                    'T_ncdm':0.71611,

                    'output':'mPk, tCl, pCl, lCl',
                    'lensing':'yes',
                    'P_k_max_h/Mpc':1.0,
                    'non_linear':'halofit'
                    })

        elif ('planck' in baseline_name and '18' in baseline_name and 'lens' in baseline_name) or 'p18l' in baseline_name.lower():
          self.set({'omega_b':2.236219e-02,
                    'omega_cdm':1.201668e-01,
                    'H0':6.726996e+01,
                    'A_s':2.102880e-09,
                    'n_s':9.661489e-01,
                    'z_reio':7.743057e+00,

                    'N_ur':2.0328,
                    'N_ncdm':1,
                    'm_ncdm':0.06,
                    'T_ncdm':0.71611,

                    'output':'mPk, tCl, pCl, lCl',
                    'lensing':'yes',
                    'P_k_max_h/Mpc':1.0,
                    'non_linear':'halofit'
                    })

        elif ('planck' in baseline_name and '18' in baseline_name) or 'p18' in baseline_name.lower():
          self.set({'omega_b':2.237064e-02,
                    'omega_cdm':1.214344e-01,
                    'H0':6.685836e+01,
                    'A_s':2.112203e-09,
                    'n_s':9.622800e-01,
                    'z_reio':7.795700e+00,

                    'N_ur':2.0328,
                    'N_ncdm':1,
                    'm_ncdm':0.06,
                    'T_ncdm':0.71611,

                    'output':'mPk, tCl, pCl, lCl',
                    'lensing':'yes',
                    'P_k_max_h/Mpc':1.0})
        else:
          raise CosmoSevereError("Unrecognized baseline case '{}'".format(baseline_name))

    def raw_cl(self, lmax=-1, nofail=False):
        """
        Unlensed CMB spectra

        Return a dictionary of the unlensed CMB spectra C_l for all modes requested in the 'output' field
        (temperature, polarisation, lensing potential, cross-correlations...)
        This function requires that the 'ouput' field contains at least on of 'tCl', 'pCl', 'lCl'.
        The returned C_l's are all dimensionless.

        If lmax is not passed (default), the spectra are returned up to the maximum multipole that was set through input parameters.
        If lmax is passed, the code first checks whether it is smaller or equal to the maximum multipole that was set through input parameters.
        If it is smaller or equal, the spectra are returned up to lmax.
        If it is larger and nofail=False (default), the function raises an error.
        If it is larger and nofail=True, the functiona asks CLASS to recompute the Cl's up to lmax, and return the result up to lmax.

        Parameters
        ----------
        lmax : int, optional
            Define the maximum l for which the C_l will be returned (inclusively).

        nofail: bool, optional
            Even if lmax is larger than expected from the user input, check and enforce the computation of the C_l in the harmonic module up to lmax.

        Returns
        -------
        cl : dict
            Dictionary that contains the unlensed power spectrum for each auto-correlation and cross-correlation spectrum,
            as well the corresponding mutipoles l. The keys of the dictionary are: 'ell' for the multipoles,
            and 'tt', 'ee', 'te', 'bb', 'pp', 'tp' for each spectrum type.
        """
        self.compute(["harmonic"])
        cdef int lmaxR

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

        # We need to be able to gracefully exit BEFORE allocating things (!)
        if not spectra:
            raise CosmoSevereError("No Cl computed")

        # We need to be able to gracefully exit BEFORE allocating things (!)
        lmaxR = self.hr.l_max_tot
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        # Now that the conditions are all checked, we can allocate and do what we want

        #temporary storage for the cls (total)
        cdef double *rcl = <double*> calloc(self.hr.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md[index_md] = <double*> calloc(self.hr.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.hr.ct_size*self.hr.ic_ic_size[index_md], sizeof(double))

        # Initialise all the needed Cls arrays
        cl = {}
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)

        success = True
        # Recover for each ell the information from CLASS
        for ell from 2<=ell<lmax+1:
            if harmonic_cl_at_l(&self.hr, ell, rcl, cl_md, cl_md_ic) == _FAILURE_:
                success = False
                break
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

        # This has to be delayed until AFTER freeing the memory
        if not success:
          raise CosmoSevereError(self.hr.error_message)

        return cl

    def lensed_cl(self, lmax=-1,nofail=False):
        """
        Lensed CMB spectra

        Return a dictionary of the lensed CMB spectra C_l for all modes requested in the 'output' field
        (temperature, polarisation, lensing potential, cross-correlations...)
        This function requires that the 'ouput' field contains at least one of 'tCl', 'pCl', 'lCl', and that 'lensing' is set to 'yes'.
        The returned spectra are all dimensionless.

        If lmax is not passed (default), the spectra are returned up to the maximum multipole that was set through input parameters.
        If lmax is passed, the code first checks whether it is smaller or equal to the maximum multipole that was set through input parameters.
        If it is smaller or equal, the spectra are returned up to lmax.
        If it is larger and nofail=False (default), the function raises an error.
        If it is larger and nofail=True, the functiona asks CLASS to recompute the Cl's up to lmax, and return the result up to lmax.

        Parameters
        ----------
        lmax : int, optional
            Define the maximum l for which the C_l will be returned (inclusively).

        nofail: bool, optional
            Even if lmax is larger than expected from the user input, check and enforce the computation of the C_l in the lensing module up to lmax.

        Returns
        -------
        cl : dict
            Dictionary that contains the lensed power spectrum for each auto-correlation and cross-correlation spectrum,
            as well the corresponding mutipoles l. The keys of the dictionary are: 'ell' for the multipoles,
            and 'tt', 'ee', 'te', 'bb', 'pp', 'tp' for each spectrum type.
        """
        self.compute(["lensing"])
        cdef int lmaxR

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

        # We need to be able to gracefully exit BEFORE allocating things (!)
        if not spectra:
            raise CosmoSevereError("No lensed Cl computed")

        # We need to be able to gracefully exit BEFORE allocating things (!)
        lmaxR = self.le.l_lensed_max
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        # Now that the conditions are all checked, we can allocate and do what we want
        cdef double *lcl = <double*> calloc(self.le.lt_size,sizeof(double))

        cl = {}
        success = True
        # Simple Cls, for temperature and polarisation, are not so big in size
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)
        for ell from 2<=ell<lmax+1:
            if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
                success = False
                break
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = lcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(lcl)

        # This has to be delayed until AFTER freeing the memory
        if not success:
          raise CosmoSevereError(self.le.error_message)

        return cl

    def density_cl(self, lmax=-1, nofail=False):
        """
        C_l spectra of density and lensing

        Return a dictionary of the spectra C_l of large scale
        structure observables, that is, density (e.g. from galaxy number count)
        and lensing potential (related to cosmic shear), for all the
        reshift bins defined by the user in input. If CMB temperature
        was also requested, this function also returns the
        density-temperature and lensing-temperature cross-correlation
        spectra.
        This function requires that the 'ouput' field contains at least one of 'dCl', 'sCl'.
        The returned spectra are all dimensionless.

        If lmax is not passed (default), the spectra are returned up to the maximum multipole that was set through input parameters.
        If lmax is passed, the code first checks whether it is smaller or equal to the maximum multipole that was set through input parameters.
        If it is smaller or equal, the spectra are returned up to lmax.
        If it is larger and nofail=False (default), the function raises an error.
        If it is larger and nofail=True, the functiona asks CLASS to recompute the Cl's up to lmax, and return the result up to lmax.

        Parameters
        ----------
        lmax : int, optional
            Define the maximum l for which the C_l will be returned (inclusively).

        nofail: bool, optional
            Even if lmax is larger than expected from the user input, check and enforce the computation of the C_l in the lensing module up to lmax.

        Returns
        -------
        cl : dict
            Dictionary that contains the power spectrum for each auto-correlation and cross-correlation spectrum,
            as well the corresponding mutipoles l. The keys of the dictionary are: 'ell' for the multipoles,
            and 'dens', 'td', 'll', 'dl', 'tl', where 'd' means density, 'l' means lensing, and 't' means temperature.
        """
        self.compute(["harmonic"])
        cdef int lmaxR

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

        # We need to be able to gracefully exit BEFORE allocating things (!)
        if not spectra:
            raise CosmoSevereError("No density Cl computed")

        # We need to be able to gracefully exit BEFORE allocating things (!)
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_lss",lmax)
                self._pars_check("output",'nCl')
                self.compute()
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        # Now that the conditions are all checked, we can allocate and do what we want
        cdef double *dcl = <double*> calloc(self.hr.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md[index_md] = <double*> calloc(self.hr.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.hr.md_size, sizeof(double*))
        for index_md in range(self.hr.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.hr.ct_size*self.hr.ic_ic_size[index_md], sizeof(double))

        cl = {}

        # For density Cls, we compute the names for each combination, which will also correspond to the size
        names = {'dd':[],'ll':[],'dl':[]}
        for index_d1 in range(self.hr.d_size):
          for index_d2 in range(index_d1, min(index_d1+self.hr.non_diag+1, self.hr.d_size)):
            names['dd'].append("dens[%d]-dens[%d]"%(index_d1+1, index_d2+1))
            names['ll'].append("lens[%d]-lens[%d]"%(index_d1+1, index_d2+1))
          for index_d2 in range(max(index_d1-self.hr.non_diag,0), min(index_d1+self.hr.non_diag+1, self.hr.d_size)):
            names['dl'].append("dens[%d]-lens[%d]"%(index_d1+1, index_d2+1))

        for elem in names:
            if elem in spectra:
                cl[elem] = {}
                for name in names[elem]:
                    cl[elem][name] = np.zeros(lmax+1, dtype=np.double)

        for elem in ['td', 'tl']:
            if elem in spectra:
                cl[elem] = np.zeros(lmax+1, dtype=np.double)

        success = True
        for ell from 2<=ell<lmax+1:
            if harmonic_cl_at_l(&self.hr, ell, dcl, cl_md, cl_md_ic) == _FAILURE_:
                success = False
                break
            if 'dd' in spectra:
                for index, name in enumerate(names['dd']):
                  cl['dd'][name][ell] = dcl[self.hr.index_ct_dd+index]
            if 'll' in spectra:
                for index, name in enumerate(names['ll']):
                  cl['ll'][name][ell] = dcl[self.hr.index_ct_ll+index]
            if 'dl' in spectra:
                for index, name in enumerate(names['dl']):
                  cl['dl'][name][ell] = dcl[self.hr.index_ct_dl+index]
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

        # This has to be delayed until AFTER freeing the memory
        if not success:
          raise CosmoSevereError(self.hr.error_message)
        return cl

    def z_of_r (self, z):
        """
        Return the comoving radius r(z) and the derivative dz/dr.

        Return the comoving radius r(z) (in units of Mpc) of an object
        seen at redshift z (dimensionless), as well and the derivative
        dz/dr (in units of 1/Mpc). The name of the function is
        misleading, it should be r_of_z. This naming inconsistency
        propagated until now and could be fixed at some point. The
        function accepts single entries or arrays of z values.

        Parameters
        ----------
        z : float array
            Redshift

        Returns
        -------
        float array
            Comoving radius

        real array
            Derivative dz/dr
        """
        self.compute(["background"])
        cdef int last_index=0 #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        r = np.zeros(len(zarr),'float64')
        dzdr = np.zeros(len(zarr),'float64')

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        i = 0
        for redshift in zarr:

            if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                free(pvecback) #manual free due to error
                raise CosmoSevereError(self.ba.error_message)

            # store r
            r[i] = pvecback[self.ba.index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[self.ba.index_bg_H]

            i += 1

        free(pvecback)

        return (r[0], dzdr[0]) if np.isscalar(z) else (r,dzdr)

    def luminosity_distance(self, z):
        """
        Return the luminosity_distance d_L(z)

        Return the luminosity distance d_L(z) in units of Mpc.
        The function accepts single entries or arrays of z values.

        Parameters
        ----------
        z : float array
            Redshift

        Returns
        -------
        float array
            Luminosity distance
        """
        self.compute(["background"])

        cdef int last_index = 0  # junk

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        lum_distance = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba, redshift, long_info,
                  inter_normal, &last_index, pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          lum_distance[iz] = pvecback[self.ba.index_bg_lum_distance]
        free(pvecback)

        return (lum_distance[0] if np.isscalar(z) else lum_distance)

    def pk(self,double k,double z):
        """
        Return the total matter power spectrum P_m(k,z)

        Return the total matter power spectrum P_m(k,z) (in Mpc**3) for a given k (in
        1/Mpc) and z. The function returns the linear power spectrum
        if the user sets 'non_linear' to 'no', and the non-linear
        power spectrum otherwise.
        This function requires that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift

        Returns
        -------
        pk : float
            Total matter power spectrum
        """
        self.compute(["fourier"])

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

    def pk_cb(self,double k,double z):
        """
        Return the CDM+baryon power spectrum P_cb(k,z)

        Return the CDM+baryon power spectrum P_cb(k,z) (in Mpc**3) for a given k (in
        1/Mpc) and z. The function returns the linear power spectrum
        if the user sets 'non_linear' to 'no', and the non-linear
        power spectrum otherwise.
        This function requires that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift

        Returns
        -------
        pk_cb : float
            CDM+baryon power spectrum
        """
        self.compute(["fourier"])

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

    def pk_lin(self,double k,double z):
        """
        Return the linear total matter power spectrum P_m(k,z)

        Return the linear total matter power spectrum P_m(k,z) (in Mpc**3) for a given k (in
        1/Mpc) and z. The function returns the linear power spectrum
        in all cases, independently of what the user sets for 'non_linear'.
        This function requires that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift

        Returns
        -------
        pk_lin : float
            Linear total matter power spectrum
        """
        self.compute(["fourier"])

        cdef double pk_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_m,&pk_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_lin

    def pk_cb_lin(self,double k,double z):
        """
        Return the linear CDM+baryon power spectrum P_cb(k,z)

        Return the linear CDM+baryon power spectrum P_cb(k,z) (in Mpc**3) for a given k (in
        1/Mpc) and z. The function returns the linear power spectrum
        in all cases, independently of what the user sets for 'non_linear'.
        This function requires that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift

        Returns
        -------
        pk_cb_lin : float
            Linear CDM+baryon power spectrum
        """
        self.compute(["fourier"])

        cdef double pk_cb_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_cb,&pk_cb_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_cb_lin

    def pk_numerical_nw(self,double k,double z):
        """
        Return the no-wiggle (smoothed) linear matter power spectrum P_m(k,z)

        Return the no-wiggle (smoothed) linear matter power spectrum P_m(k,z) (in Mpc**3) for a given k (in
        1/Mpc) and z. A smoothing algorithm infers this spectrum from the full linear matter power spectrum.
        This function requires that the 'ouput' field contains at least 'mPk', and 'numerical_nowiggle' is set to 'yes'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift

        Returns
        -------
        pk_numerical_nw : float
            No-wiggle linear matter power spectrum
        """
        self.compute(["fourier"])

        cdef double pk_numerical_nw

        if (self.fo.has_pk_numerical_nowiggle == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must set `numerical_nowiggle` to `yes` in input")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_numerical_nowiggle,k,z,0,&pk_numerical_nw,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_numerical_nw

    def pk_analytic_nw(self,double k):
        """
        Return an analytic approximation to the no-wiggle linear matter power spectrum P_m(k) at z=0

        Return an analytic approximation to the linear matter power spectrum P_m(k) (in Mpc**3) for a given k (in
        1/Mpc) and at z=0, without BAO features. The calculation is based on the Eisenstein & Hu fitting formulas.
        This function requires that the 'ouput' field contains at least 'mPk', and 'analytic_nowiggle' is set to 'yes'.

        Parameters
        ----------
        k : float
            Wavenumber

        Returns
        -------
        pk_analytic_nw : float
            Analytic approximation to no-wiggle linear matter power spectrum
        """
        self.compute(["fourier"])

        cdef double pk_analytic_nw

        if (self.fo.has_pk_analytic_nowiggle == _FALSE_):
            raise CosmoSevereError("No analytic nowiggle spectrum computed. You must set `analytic_nowiggle` to `yes` in input")

        if fourier_pk_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_analytic_nowiggle,k,0.,self.fo.index_pk_m,&pk_analytic_nw,NULL)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_analytic_nw

    def get_pk(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """
        Return the total matter power spectrum P_m(k,z) for a 3D array of k and a 1D array of z

        Return the total matter power spectrum P_m(k,z) (in Mpc**3)
        for a 3D array of k values (in 1/Mpc) and z values. The array
        of k values must be indexed as k[index_k,index_z,index_mu],
        the one of z values as z[index_z]. This function is useful in
        the context of likelihoods describing spectroscopic galaxy
        redshift, when there is a different k-list for each redhsift
        and each value of the angle mu between the line-of-sight and
        the measured two-point correlation function. The function
        returns a grid of values for the linear power spectrum if the
        user sets 'non_linear' to 'no', and for the non-linear power
        spectrum otherwise.  This function requires that the 'ouput'
        field contains at least 'mPk'.

        Parameters
        ----------
        k : float array
            Wavenumber array indexed as k[index_k,index_z,index_mu]
        z : float array
            Redshift array indexed as z[index_z]
        k_size : int
            Number of k values for each index_z and index_mu
        z_size : int
            Number of redshift values
        mu_size : int
            Number of k values for each index_k and index_z

        Returns
        -------
        pk : float
            Grid of total matter power spectrum indexed as pk[index_k,index_z,index_mu]
        """
        self.compute(["fourier"])

        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in range(k_size):
            for index_z in range(z_size):
                for index_mu in range(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """
        Return the CDM+baryon power spectrum P_cb(k,z) for a 3D array of k and a 1D array of z

        Return the CDM+baryon power spectrum P_m(k,z) (in Mpc**3)
        for a 3D array of k values (in 1/Mpc) and z values. The array
        of k values must be indexed as k[index_k,index_z,index_mu],
        the one of z values as z[index_z]. This function is useful in
        the context of likelihoods describing spectroscopic galaxy
        redshift, when there is a different k-list for each redhsift
        and each value of the angle mu between the line-of-sight and
        the measured two-point correlation function. The function
        returns a grid of values for the linear power spectrum if the
        user sets 'non_linear' to 'no', and for the non-linear power
        spectrum otherwise.  This function requires that the 'ouput'
        field contains at least 'mPk'.

        Parameters
        ----------
        k : float array
            Wavenumber array indexed as k[index_k,index_z,index_mu]
        z : float array
            Redshift array indexed as z[index_z]
        k_size : int
            Number of k values for each index_z and index_mu
        z_size : int
            Number of redshift values
        mu_size : int
            Number of k values for each index_k and index_z

        Returns
        -------
        pk_cb : float
            Grid of CDM+baryon power spectrum indexed as pk_cb[index_k,index_z,index_mu]
        """
        self.compute(["fourier"])

        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in range(k_size):
            for index_z in range(z_size):
                for index_mu in range(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    def get_pk_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """
        Return the linear total matter power spectrum P_m(k,z) for a 3D array of k and a 1D array of z

        Return the linear total matter power spectrum P_m(k,z) (in
        Mpc**3) for a 3D array of k values (in 1/Mpc) and z
        values. The array of k values must be indexed as
        k[index_k,index_z,index_mu], the one of z values as
        z[index_z]. This function is useful in the context of
        likelihoods describing spectroscopic galaxy redshift, when
        there is a different k-list for each redhsift and each value
        of the angle mu between the line-of-sight and the measured
        two-point correlation function. The function always returns a
        grid of values for the linear power spectrum, independently of
        what the user sets for 'non_linear'.  This function requires
        that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float array
            Wavenumber array indexed as k[index_k,index_z,index_mu]
        z : float array
            Redshift array indexed as z[index_z]
        k_size : int
            Number of k values for each index_z and index_mu
        z_size : int
            Number of redshift values
        mu_size : int
            Number of k values for each index_k and index_z

        Returns
        -------
        pk : float
            Grid of linear total matter power spectrum indexed as pk[index_k,index_z,index_mu]
        """
        self.compute(["fourier"])

        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in range(k_size):
            for index_z in range(z_size):
                for index_mu in range(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """
        Return the linear CDM+baryon power spectrum P_cb(k,z) for a 3D array of k and a 1D array of z

        Return the linear CDM+baryon matter power spectrum P_m(k,z) (in
        Mpc**3) for a 3D array of k values (in 1/Mpc) and z
        values. The array of k values must be indexed as
        k[index_k,index_z,index_mu], the one of z values as
        z[index_z]. This function is useful in the context of
        likelihoods describing spectroscopic galaxy redshift, when
        there is a different k-list for each redhsift and each value
        of the angle mu between the line-of-sight and the measured
        two-point correlation function. The function always returns a
        grid of values for the linear power spectrum, independently of
        what the user sets for 'non_linear'.  This function requires
        that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float array
            Wavenumber array indexed as k[index_k,index_z,index_mu]
        z : float array
            Redshift array indexed as z[index_z]
        k_size : int
            Number of k values for each index_z and index_mu
        z_size : int
            Number of redshift values
        mu_size : int
            Number of k values for each index_k and index_z

        Returns
        -------
        pk_cb : float
            Grid of linear CDM+baryon power spectrum indexed as pk_cb[index_k,index_z,index_mu]
        """
        self.compute(["fourier"])

        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in range(k_size):
            for index_z in range(z_size):
                for index_mu in range(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    def get_pk_all(self, k, z, nonlinear = True, cdmbar = False, z_axis_in_k_arr = 0, interpolation_kind='cubic'):
        """
        Return different versions of the power spectrum P(k,z) for arrays of k and z with arbitrary shapes

        Return different versions of the power spectrum P_m(k,z) (in Mpc**3)
        for arrays of k values (in 1/Mpc) and z values with arbitrary shape.
        The optional argument nonlinear can be set to True for non-linear (default) or False for non-linear.
        The optional argument cdmbar can be set to False for total matter (default) or True for CDM+baryon. It makes a difference only in presence of massive neutrinos.
        For multi-dimensional k-arrays, the function assumes that one of the dimensions is the z-axis.
        The optional argument z_axis_in_k_arr specifies the integer position of the z_axis within the n-dimensional k array (default:0).
        This function requires that the 'ouput' field contains at least 'mPk'.
        The function returns a grid of values for the power spectrum, with the same shape as the input k grid.

        Parameters
        ----------
        k : float array
            Wavenumber array, arbitrary shape
        z : float array
            Redshift array, arbitrary shape
        nonlinear: bool, optional
            Whether to return the non-linear spectrum instead of the linear one
        cdmbar : bool, optional
            Whether to return the CDM+baryon spectrum instead of the total matter one
        z_axis_in_k_arr : int, optional
            Position of the z_axis within the k array
        interpolation_kind : str, optional
            Flag for the interpolation kind, to be understood by interp1d() function of scipy.interpolate module

        Returns
        -------
        out_pk : float
            Grid of power spectrum with the same shape as input k grid
        """
        self.compute(["fourier"])

        # z_axis_in_k_arr specifies the integer position of the z_axis wihtin the n-dimensional k_arr
        # Example: 1-d k_array -> z_axis_in_k_arr = 0
        # Example: 3-d k_array with z_axis being the first axis -> z_axis_in_k_arr = 0
        # Example: 3-d k_array with z_axis being the last axis  -> z_axis_in_k_arr = 2

        # 1) Define some utilities
        # Is the user asking for a valid cdmbar?
        ispkcb = cdmbar and not (self.ba.Omega0_ncdm_tot == 0.)

        # Allocate the temporary k/pk array used during the interaction with the underlying C code
        cdef np.float64_t[::1] pk_out = np.empty(self.fo.k_size, dtype='float64')
        k_out = np.asarray(<np.float64_t[:self.fo.k_size]> self.fo.k)

        # Define a function that can write the P(k) for a given z into the pk_out array
        def _write_pk(z,islinear,ispkcb):
          if fourier_pk_at_z(&self.ba,&self.fo,linear,(pk_linear if islinear else pk_nonlinear),z,(self.fo.index_pk_cb if ispkcb else self.fo.index_pk_m),&pk_out[0],NULL)==_FAILURE_:
              raise CosmoSevereError(self.fo.error_message)

        # Check what kind of non-linear redshift there is
        if nonlinear:
          if self.fo.index_tau_min_nl == 0:
            z_max_nonlinear = np.inf
          else:
            z_max_nonlinear = self.z_of_tau(self.fo.tau[self.fo.index_tau_min_nl])
        else:
          z_max_nonlinear = -1.

        # Only get the nonlinear function where the nonlinear treatment is possible
        def _islinear(z):
          if z > z_max_nonlinear or (self.fo.method == nl_none):
            return True
          else:
            return False

        # A simple wrapper for writing the P(k) in the given location and interpolating it
        def _interpolate_pk_at_z(karr,z):
          _write_pk(z,_islinear(z),ispkcb)
          interp_func = interp1d(k_out,np.log(pk_out),kind=interpolation_kind,copy=True)
          return np.exp(interp_func(karr))

        # 2) Check if z array, or z value
        if not isinstance(z,(list,np.ndarray)):
            # Only single z value was passed -> k could still be an array of arbitrary dimension
            if not isinstance(k,(list,np.ndarray)):
                # Only single z value AND only single k value -> just return a value
                # This iterates over ALL remaining dimensions
                return ((self.pk_cb if ispkcb else self.pk) if not _islinear(z) else (self.pk_cb_lin if ispkcb else self.pk_lin))(k,z)
            else:
                k_arr = np.array(k)
                result = _interpolate_pk_at_z(k_arr,z)
                return result

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
                out_pk[index_z] = _interpolate_pk_at_z(k_arr[index_z],z[index_z])
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
                    out_pk_at_z = _interpolate_pk_at_z(k_arr_at_z,z[index_z])
                    out_pk.append(out_pk_at_z)
                return out_pk

            # 3.3) If there is a single-dimensional k-array
            # The user passed a z-array, but only a 1-d k array
            # Assume thus, that the k array should be reproduced for all z
            out_pk = np.empty((len(z_arr),len(k_arr)))
            for index_z in range(len(z_arr)):
                out_pk[index_z] = _interpolate_pk_at_z(k_arr,z_arr[index_z])
            return out_pk

    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species = False, h_units=False):
        """
        Return different versions of the power spectrum P(k,z), together with the corresponding list of k and z values

        Return different versions of the power spectrum P(k,z) (in Mpc**3) together with the corresponding list
        k values (in 1/Mpc) and z values. The values are those defined internally in CLASS before interpolation. This function is useful for building interpolators.
        The optional argument nonlinear can be set to True for non-linear (default) or False for non-linear.
        The optional argument only_clustering_species can be set to False for total matter (default) or True for CDM+baryon. It makes a difference only in presence of massive neutrinos.
        The optional argument h_units, set to False by default, can be set to True to have P(k,z) in (Mpc/h)**3 and k in (h/Mpc).
        This function requires that the 'ouput' field contains at least 'mPk'.
        The function returns a 2D grid of values for the power spectrum.

        Parameters
        ----------
        nonlinear : bool, optional
            Whether to return the non-linear spectrum instead of the linear one
        only_clustering_species : bool, optional
            Whether to return the CDM+baryon spectrum instead of the total matter one
        h_units : bool, optional
            Whether to use Mpc/h and h/Mpc units instead of Mpc and 1/Mpc

        Returns
        -------
        pk : float array
            Grid of power spectrum indexed as pk[index_k, index_z]
        k : float array
            Vector of wavenumbers indexed as k[index_k]
        z : float array
            Vector of redshifts indexed as k[index_z]
        """
        self.compute(["fourier"])

        cdef np.ndarray[DTYPE_t,ndim=2] pk = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.fo.k_size_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.fo.ln_tau_size),'float64')
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
        # For nonlinear, we have to additionally cut out the linear values

        if self.fo.ln_tau_size == 1:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in range(self.fo.ln_tau_size):
                if index_tau == self.fo.ln_tau_size-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.fo.ln_tau[index_tau]))

        # check consitency of the list of redshifts

        if nonlinear == True:
            # Check highest value of z at which nl corrections could be computed.
            # In the table tau_sampling it corresponds to index: self.fo.index_tau_min_nl
            z_max_nonlinear = self.z_of_tau(self.fo.tau[self.fo.index_tau_min_nl])

            # Check highest value of z in the requested output.
            z_max_requested = z[0]

            # The first z must be larger or equal to the second one, that is,
            # the first index must be smaller or equal to the second one.
            # If not, raise and error.
            if (z_max_requested > z_max_nonlinear and self.fo.index_tau_min_nl>0):
                raise CosmoSevereError("get_pk_and_k_and_z() is trying to return P(k,z) up to z_max=%e (the redshift range of computed pk); but the input parameters sent to CLASS (in particular ppr->nonlinear_min_k_max=%e) were such that the non-linear P(k,z) could only be consistently computed up to z=%e; increase the precision parameter 'nonlinear_min_k_max', or only obtain the linear pk"%(z_max_requested,self.pr.nonlinear_min_k_max,z_max_nonlinear))

        # get list of k

        if h_units:
            units=1./self.ba.h
        else:
            units=1

        for index_k in range(self.fo.k_size_pk):
            k[index_k] = self.fo.k[index_k]*units

        # get P(k,z) array

        for index_tau in range(self.fo.ln_tau_size):
            for index_k in range(self.fo.k_size_pk):
                if nonlinear == True:
                    pk[index_k, index_tau] = np.exp(self.fo.ln_pk_nl[index_pk][index_tau * self.fo.k_size + index_k])
                else:
                    pk[index_k, index_tau] = np.exp(self.fo.ln_pk_l[index_pk][index_tau * self.fo.k_size + index_k])

        return pk, k, z

    def get_transfer_and_k_and_z(self, output_format='class', h_units=False):
        """
        Return a grid of transfer function of various perturbations T_i(k,z) arranged in a dictionary, together with the corresponding list of k and z values

        Return a grid of transfer function of various density and/or velocity perturbations T_i(k,z) arranged in a dictionary, together with the corresponding list of k values (in 1/Mpc) and z values. The values are those defined internally in CLASS before interpolation. This function is useful for building interpolators.
        The optional argument output_format can be set to 'class' (default) or 'camb'. With output_format='class', all transfer functions will be normalised to 'curvature R=1' at initial time and are dimensionless. With output_format='camb', they are normalised to 'curvature R = -1/k^2' like in CAMB, and they have units of Mpc**2. Then, 'dTk' must be in the input: the CAMB format only outputs density transfer functions.
        When sticking to output_format='class', you also get the newtonian metric fluctuations phi and psi.
        If you set the CLASS input parameter 'extra_metric_transfer_functions' to 'yes',
        you get additional metric fluctuations in the synchronous and N-body gauges.
        The optional argument h_units, set to False by default, can be set to True to have k in (h/Mpc).
        This function requires that the 'ouput' field contains at least one of 'dTk' (for density transfer functions) or 'vTk' (for velocity transfer functions). With output_format='camb', 'dTk' must be in the input: the CAMB format only outputs density transfer functions.
        The function returns a dictionary with, for each species, a 2D grid of transfer functions.

        Parameters
        ----------
        output_format : str, optional
            Format transfer functions according to 'class' or 'camb' convention
        h_units : bool, optional
            Whether to use h/Mpc units for k instead of 1/Mpc

        Returns
        -------
        tk : dict
            Dictionary containing all transfer functions.
            For instance, the grid of values of 'd_c' (= delta_cdm) is available in tk['d_c']
            All these grids are indexed as [index_k,index,z], for instance tk['d_c'][index_k,index,z]
        k : float array
            Vector of wavenumbers indexed as k[index_k]
        z : float array
            Vector of redshifts indexed as k[index_z]
        """
        self.compute(["transfer"])

        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.pt.k_size_pk),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.pt.ln_tau_size),'float64')
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
            free(titles) # manual free due to error
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]

        free(titles)

        number_of_titles = len(names)

        # get list of redshifts
        # the ln(times) of interest are stored in self.fo.ln_tau[index_tau]

        if self.pt.ln_tau_size == 1:
            raise CosmoSevereError("You ask classy to return an array of T_x(k,z) values, but the input parameters sent to CLASS did not require any transfer function calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in range(self.pt.ln_tau_size):
                if index_tau == self.pt.ln_tau_size-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.pt.ln_tau[index_tau]))

        # get list of k

        if h_units:
            units=1./self.ba.h
        else:
            units=1

        k_size = self.pt.k_size_pk
        for index_k in range(k_size):
            k[index_k] = self.pt.k[index_md][index_k]*units

        # create output dictionary

        tk = {}
        for index_type,name in enumerate(names):
            if index_type > 0:
                tk[name] = np.zeros((k_size, len(z)),'float64')

        # allocate the vector in wich the transfer functions will be stored temporarily for all k and types at a given z
        data = <double*>malloc(sizeof(double)*number_of_titles*self.pt.k_size[index_md])

        # get T(k,z) array

        for index_tau in range(len(z)):
            if perturbations_output_data_at_index_tau(&self.ba, &self.pt, outf, index_tau, number_of_titles, data)==_FAILURE_:
                free(data) # manual free due to error
                raise CosmoSevereError(self.pt.error_message)

            for index_type,name in enumerate(names):
                if index_type > 0:
                    for index_k in range(k_size):
                        tk[name][index_k, index_tau] = data[index_k*number_of_titles+index_type]

        free(data)
        return tk, k, z

    def get_Weyl_pk_and_k_and_z(self, nonlinear=False, h_units=False):
        """
        Return the power spectrum of the perturbation [k^2*(phi+psi)/2](k,z), together with the corresponding list of k and z values

        Return the power spectrum of the perturbation [k^2*(phi+psi)/2](k,z), where (phi+psi)/2 is the Weyl potential transfer function, together with the corresponding list of k values (in 1/Mpc) and z values. The values are those defined internally in CLASS before interpolation. This function is useful for building interpolators.
        Note that this function first gets P(k,z) and then corrects the output
        by the ratio of transfer functions [k^2(phi+psi)/d_m]^2.
        The optional argument nonlinear can be set to True for non-linear (default) or False for non-linear.
        The optional argument h_units, set to False by default, can be set to True to have k in (h/Mpc).
        This function requires that the 'ouput' field contains at least one of 'mTk' or 'vTk'.
        The function returns a 2D grid of values for the Weyl potential.

        Parameters
        ----------
        nonlinear : bool, optional
            Whether to return the non-linear spectrum instead of the linear one
        h_units : bool, optional
            Whether to use Mpc/h and h/Mpc units instead of Mpc and 1/Mpc

        Returns
        -------
        Weyl_pk : float array
            Grid of Weyl_pk indexed as Weyl_pk[index_k, index_z]
        k : float array
            Vector of wavenumbers indexed as k[index_k]
        z : float array
            Vector of redshifts indexed as k[index_z]
        """
        self.compute(["fourier"])

        cdef np.ndarray[DTYPE_t,ndim=2] pk = np.zeros((self.fo.k_size_pk,self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] k4 = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] phi = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] psi = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] d_m = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=2] Weyl_pk = np.zeros((self.fo.k_size_pk, self.fo.ln_tau_size),'float64')

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
        for index_z in range(self.fo.ln_tau_size):
            k4[:,index_z] = k**4

        # rescale total matter power spectrum to get the Weyl power spectrum times k**4
        # (the latter factor is just a convention. Since there is a factor k**2 in the Poisson equation,
        # this rescaled Weyl spectrum has a shape similar to the matter power spectrum).
        Weyl_pk = pk * ((phi+psi)/2./d_m)**2 * k4

        return Weyl_pk, k, z

    def sigma(self,R,z, h_units = False):
        """
        Return sigma (total matter) for radius R and redhsift z

        Return sigma, the root mean square (rms) of the relative density fluctuation of
        total matter (dimensionless) in spheres of radius R at redshift z. If R= 8 h/Mpc this will be the usual sigma8(z).
        If the user passes a single value of R and z, the function returns sigma(R,z).
        If the use rpasses an array of R and/or z values, the function sets the shape of the returned grid of sigma(R,z) accordingly.

        If h_units = False (default), R is in unit of Mpc and sigma8 is obtained for R = 8/h.
        If h_units = True, R is in unit of Mpc/h and sigma8 is obtained for R = 8.

        This function requires that the 'ouput' field contains at
        least 'mPk' and that 'P_k_max_h/Mpc' or 'P_k_max_1/Mpc' are
        such that k_max is bigger or equal to 1 h/Mpc.

        Parameters
        ----------
        R : float
            Radius of the spheres in which the rms is computed (single value or array)
        z : float
            Redshift (single value or array)
        h_units : bool, optional
            Whether to use Mpc/h instead of Mpc for R

        Returns
        -------
        sigma : float
            Rms of density fluctuation of total matter (single value or grid of values)
        """
        self.compute(["fourier"])

        cdef double sigma

        zarr = np.atleast_1d(z).astype(np.float64)
        Rarr = np.atleast_1d(R).astype(np.float64)

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        R_in_Mpc = (Rarr if not h_units else Rarr/self.ba.h)

        pairs = np.array(np.meshgrid(zarr,R_in_Mpc)).T.reshape(-1,2)

        sigmas = np.empty(pairs.shape[0])
        for ip, pair in enumerate(pairs):
          if fourier_sigmas_at_z(&self.pr,&self.ba,&self.fo,pair[1],pair[0],self.fo.index_pk_m,out_sigma,&sigma)==_FAILURE_:
              raise CosmoSevereError(self.fo.error_message)
          sigmas[ip] = sigma

        return (sigmas[0] if (np.isscalar(z) and np.isscalar(R)) else np.squeeze(sigmas.reshape(len(zarr),len(Rarr))))

    def sigma_cb(self,double R,double z, h_units = False):
        """
        Return sigma_cb (CDM+baryon) for radius R and redhsift z

        Return sigma_cb, the root mean square (rms) of the relative density fluctuation of
        CDM+baryon (dimensionless) in spheres of radius R at redshift z. If R= 8 h/Mpc this will be the usual sigma8_cb(z).
        If the user passes a single value of R and z, the function returns sigma_cb(R,z).
        If the user passes an array of R and/or z values, the function sets the shape of the returned grid of sigma_cb(R,z) accordingly.

        If h_units = False (default), R is in unit of Mpc and sigma8_cb is obtained for R = 8/h.
        If h_units = True, R is in unit of Mpc/h and sigma8_cb is obtained for R = 8.

        This function requires that the 'ouput' field contains at
        least 'mPk' and that 'P_k_max_h/Mpc' or 'P_k_max_1/Mpc' are
        such that k_max is bigger or equal to 1 h/Mpc.

        Parameters
        ----------
        R : float
            Radius of the spheres in which the rms is computed (single value or array)
        z : float
            Redshift (single value or array)
        h_units : bool, optional
            Whether to use Mpc/h instead of Mpc for R

        Returns
        -------
        sigma_cb : float
            Rms of density fluctuation of CDM+baryon (single value or grid of values)
        """
        self.compute(["fourier"])

        cdef double sigma_cb

        zarr = np.atleast_1d(z).astype(np.float64)
        Rarr = np.atleast_1d(R).astype(np.float64)

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.fo.has_pk_cb == _FALSE_):
            raise CosmoSevereError("sigma_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        R_in_Mpc = (Rarr if not h_units else Rarr/self.ba.h)

        pairs = np.array(np.meshgrid(zarr,R_in_Mpc)).T.reshape(-1,2)

        sigmas_cb = np.empty(pairs.shape[0])
        for ip, pair in enumerate(pairs):
          if fourier_sigmas_at_z(&self.pr,&self.ba,&self.fo,R,z,self.fo.index_pk_cb,out_sigma,&sigma_cb)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)
          sigmas_cb[ip] = sigma_cb

        return (sigmas_cb[0] if (np.isscalar(z) and np.isscalar(R)) else np.squeeze(sigmas_cb.reshape(len(zarr),len(Rarr))))

    def pk_tilt(self,double k,double z):
        """
        Return the logarithmic slope of the linear matter power spectrum at k and z

        Return the logarithmic slope of the linear matter power
        spectrum, d ln P / d ln k (dimensionless), at a given
        wavenumber k (units of 1/Mpc) and redshift z.

        This function requires that the 'ouput' field contains at
        least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift

        Returns
        -------
        pk_tilt : float
            Logarithmic slope
        """
        self.compute(["fourier"])

        cdef double pk_tilt

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get pk_tilt(k,z) you must add mPk to the list of outputs.")

        if (k < self.fo.k[1] or k > self.fo.k[self.fo.k_size-2]):
            raise CosmoSevereError("In order to get pk_tilt at k=%e 1/Mpc, you should compute P(k,z) in a wider range of k's"%k)

        if fourier_pk_tilt_at_k_and_z(&self.ba,&self.pm,&self.fo,pk_linear,k,z,self.fo.index_pk_total,&pk_tilt)==_FAILURE_:
            raise CosmoSevereError(self.fo.error_message)

        return pk_tilt

    def angular_distance(self, z):
        """
        Return the angular diameter distance at z

        Return the angular diameter distance (units of Mpc) at redshift z.
        If the user passes a single value of z, the function returns angular_distance(z).
        If the user passes an array of z values, the function sets the shape of angular_distance(z) accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        angular_distance : float
            Angular distance (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        D_A = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #Manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          D_A[iz] = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return (D_A[0] if np.isscalar(z) else D_A)

    def angular_distance_from_to(self, z1, z2):
        """
        Return the angular diameter distance of object at z2 as seen by observer at z1

        Return the angular diameter distance (units of Mpc) of object
        at z2 as seen by observer at z1, that is,
        sin_K((chi2-chi1)*np.sqrt(|k|))/np.sqrt(|k|)/(1+z2).  If
        z1>z2, returns zero.

        Parameters
        ----------
        z1 : float
            Observer redshift
        z2 : float
            Source redshift

        Returns
        -------
        angular_distance_from_to : float
            Angular distance
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        if z1>=z2:
            return 0.

        else:
            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            if background_at_z(&self.ba,z1,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                free(pvecback) #manual free due to error
                raise CosmoSevereError(self.ba.error_message)

            # This is the comoving distance to object at z1
            chi1 = pvecback[self.ba.index_bg_conf_distance]

            if background_at_z(&self.ba,z2,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                free(pvecback) #manual free due to error
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
        Return the comoving distance at z

        Return the comoving distance (units of Mpc) at redshift z.
        If the user passes a single value of z, the function returns comoving_distance(z).
        If the user passes an array of z values, the function sets the shape of comoving_distance(z) accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        r : float
            Comoving distance (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        r = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          r[iz] = pvecback[self.ba.index_bg_conf_distance]

        free(pvecback)

        return (r[0] if np.isscalar(z) else r)

    def scale_independent_growth_factor(self, z):
        """
        Return the scale-independent growth factor at z

        Return the scale-independent growth factor D(z) for CDM
        perturbations (dimensionless).  CLASS finds this quantity by
        solving the background equation D'' + a H D' + 3/2 a^2 \rho_M
        D = 0, where ' stands for a derivative w.r.t conformal time. D
        is normalized to 1 today. By default, in this equation, rho_M
        is obtained by summing over baryons and cold dark matter
        (possibly including interacting CDM), but not over
        rho_ncdm. This means that in models with massive neutrinos D
        is a approximation for the growth factor in the small-scale
        (large-k) limit (the limit in which neutrinos free stream).

        Parameters
        ----------
        z : float
            Desired redshift

        Returns
        -------
        D : float
            Scale-independent growth factor D(z)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        D = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          D[iz] = pvecback[self.ba.index_bg_D]

        free(pvecback)

        return (D[0] if np.isscalar(z) else D)

    def scale_independent_growth_factor_f(self, z):
        """
        Return the scale-independent growth rate at z

        Return the scale-independent growth rate f(z) = (d ln D / d ln
        a) for CDM perturbations (dimensionless). CLASS finds D by
        solving the background equation D'' + a H D' + 3/2 a^2 \rho_M
        D = 0, where ' stands for a derivative w.r.t conformal
        time. By default, in this equation, rho_M is obtained by
        summing over baryons and cold dark matter (possibly including
        interacting CDM), but not over rho_ncdm. This means that in
        models with massive neutrinos f is a approximation for the
        growth rate in the small-scale (large-k) limit (the limit in
        which neutrinos free stream).

        Parameters
        ----------
        z : float
            Desired redshift

        Returns
        -------
        f : float
            Scale-independent growth rate f(z)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        f = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          f[iz] = pvecback[self.ba.index_bg_f]

        free(pvecback)

        return (f[0] if np.isscalar(z) else f)

    def scale_dependent_growth_factor_f(self, k, z, h_units=False, nonlinear=False):
        """
        Return the scale-dependent growth rate of total matter at z

        Return the scale-dependent growth rate for total matter perturbations (dimensionless):

        f(k,z) = (d ln delta_m / d ln a)

        The growth rate is actually inferred from the matter power
        spectrum rather than from the delta_m perturbation, as:

        f(k,z)= 1/2 * [d ln P_m(k,a) / d ln a] = - 0.5 * (1+z) * [d ln P_m(k,z) / d z]

        where P_m(k,z) is the total matter power spectrum. This is
        evaluated using a numerical derivative computed with the
        function UnivariateSpline of the module scipy.interpolate
        (with option s=0 in UnivariateSpline).

        If h_units = False (default), k is in unit of 1/Mpc.
        If h_units = True, k is in unit of h/Mpc.

        If nonlinear=False (default), estimate f from the linear power spectrum.
        If nonlinear=True, estimate f from the non-linear power spectrum.

        This function requires that the 'ouput' field contains at
        least 'mPk'. It also requires that the user asks for P(k,z) at redshifts
        larger than 0, or sets 'z_max_pk' to a non-zero value.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift
        h_units : bool, optional
            Whether to pass k in units of h/Mpc instead of 1/Mpc
        nonlinear : bool, optional
            Whether to compute f from the non-linear power spectrum instead of the linear one

        Returns
        -------
        f : float
            Scale-dependent growth rate f(k,z)
        """
        self.compute(["fourier"])

        # build array of z values at wich P(k,z) was pre-computed by class (for numerical derivative)
        # check that P(k,z) was stored at different zs
        if self.fo.ln_tau_size > 1:
            # check that input z is in stored range
            z_max = self.z_of_tau(np.exp(self.fo.ln_tau[0]))
            if (z<0) or (z>z_max):
                raise CosmoSevereError("You asked for f(k,z) at a redshift %e outside of the computed range [0,%e]"%(z,z_max))
            # create array of zs in growing z order (decreasing tau order)
            z_array = np.empty(self.fo.ln_tau_size)
            # first redshift is exactly zero
            z_array[0]=0.
            # next values can be inferred from ln_tau table
            if (self.fo.ln_tau_size>1):
                for i in range(1,self.fo.ln_tau_size):
                    z_array[i] = self.z_of_tau(np.exp(self.fo.ln_tau[self.fo.ln_tau_size-1-i]))
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

    def scale_dependent_growth_factor_f_cb(self, k, z, h_units=False, nonlinear=False):
        """
        Return the scale-dependent growth rate of CDM+baryon at z

        Return the scale-dependent growth rate for CDM+baryon perturbations (dimensionless):

        f(k,z) = (d ln delta_cb / d ln a)

        The growth rate is actually inferred from the matter power
        spectrum rather than from the delta_cb perturbation, as:

        f(k,z)= 1/2 * [d ln P_cb(k,a) / d ln a] = - 0.5 * (1+z) * [d ln P_cb(k,z) / d z]

        where P_cb(k,z) is the CDM+baryon power spectrum. This is
        evaluated using a numerical derivative computed with the
        function UnivariateSpline of the module scipy.interpolate
        (with option s=0 in UnivariateSpline).

        If h_units = False (default), k is in unit of 1/Mpc.
        If h_units = True, k is in unit of h/Mpc.

        If nonlinear=False (default), estimate f from the linear power spectrum.
        If nonlinear=True, estimate f from the non-linear power spectrum.

        This function requires that the 'ouput' field contains at
        least 'mPk'. It also requires that the user asks for P(k,z) at redshifts
        larger than 0, or sets 'z_max_pk' to a non-zero value.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift
        h_units : bool, optional
            Whether to pass k in units of h/Mpc instead of 1/Mpc
        nonlinear : bool, optional
            Whether to compute f from the non-linear power spectrum instead of the linear one

        Returns
        -------
        f : float
            Scale-dependent growth rate f_cb(k,z)
        """
        self.compute(["fourier"])

        # build array of z values at wich P_cb(k,z) was pre-computed by class (for numerical derivative)
        # check that P_cb(k,z) was stored at different zs
        if self.fo.ln_tau_size > 1:
            # check that input z is in stored range
            z_max = self.z_of_tau(np.exp(self.fo.ln_tau[0]))
            if (z<0) or (z>z_max):
                raise CosmoSevereError("You asked for f_cb(k,z) at a redshift %e outside of the computed range [0,%e]"%(z,z_max))
            # create array of zs in growing z order (decreasing tau order)
            z_array = np.empty(self.fo.ln_tau_size)
            # first redshift is exactly zero
            z_array[0]=0.
            # next values can be inferred from ln_tau table
            if (self.fo.ln_tau_size>1):
                for i in range(1,self.fo.ln_tau_size):
                    z_array[i] = self.z_of_tau(np.exp(self.fo.ln_tau[self.fo.ln_tau_size-1-i]))
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
            Pk = self.pk_cb(k,z)
            for iz, zval in enumerate(z_array):
                Pk_array[iz] = self.pk_cb(k,zval)
        else:
            Pk = self.pk_cb_lin(k,z)
            for iz, zval in enumerate(z_array):
                Pk_array[iz] = self.pk_cb_lin(k,zval)

        # Compute derivative (d ln P / d ln z)
        dPkdz = UnivariateSpline(z_array,Pk_array,s=0).derivative()(z)

        # Compute growth factor f
        f = -0.5*(1+z)*dPkdz/Pk

        return f

    def scale_dependent_growth_factor_D(self, k, z, h_units=False, nonlinear=False):
        """
        Return the scale-dependent growth factor of total matter at z

        Return the scale-dependent growth factor D for total matter perturbations (dimensionless).

        The growth factor is actually inferred from the matter power
        spectrum rather than from the delta_m perturbation, as:

        D(k,z) = sqrt[ P_m(k,z)/P_m(k,0) ]

        where P_m(k,z) is the total matter power spectrum.

        If h_units = False (default), k is in unit of 1/Mpc.
        If h_units = True, k is in unit of h/Mpc.

        If nonlinear=False (default), estimate f from the linear power spectrum.
        If nonlinear=True, estimate f from the non-linear power spectrum.

        This function requires that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift
        h_units : bool, optional
            Whether to pass k in units of h/Mpc instead of 1/Mpc
        nonlinear : bool, optional
            Whether to compute f from the non-linear power spectrum instead of the linear one

        Returns
        -------
        D : float
            Scale-dependent growth factor D(k,z)
        """
        self.compute(["fourier"])

        # if needed, convert k to units of 1/Mpc
        if h_units:
            k = k*self.ba.h

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

        # Get P(k,z) and P(k,0)
        if use_pk_lin == False:
            Pk = self.pk(k,z)
            Pk0 = self.pk(k,0)
        else:
            Pk = self.pk_lin(k,z)
            Pk0 = self.pk_lin(k,0)

        # Compute growth factor f
        D = np.sqrt(Pk/Pk0)

        return D

    def scale_dependent_growth_factor_D_cb(self, k, z, h_units=False, nonlinear=False):
        """
        Return the scale-dependent growth factor of CDM+baryon at z

        Return the scale-dependent growth factor D_cb for CDM+baryon perturbations (dimensionless).

        The growth factor is actually inferred from the matter power
        spectrum rather than from the delta_m perturbation, as:

        D_cb(k,z) = sqrt[ P_cb(k,z)/P_cb(k,0) ]

        where P_cb(k,z) is the total matter power spectrum.

        If h_units = False (default), k is in unit of 1/Mpc.
        If h_units = True, k is in unit of h/Mpc.

        If nonlinear=False (default), estimate f from the linear power spectrum.
        If nonlinear=True, estimate f from the non-linear power spectrum.

        This function requires that the 'ouput' field contains at least 'mPk'.

        Parameters
        ----------
        k : float
            Wavenumber
        z : float
            Redshift
        h_units : bool, optional
            Whether to pass k in units of h/Mpc instead of 1/Mpc
        nonlinear : bool, optional
            Whether to compute f from the non-linear power spectrum instead of the linear one

        Returns
        -------
        D_cb : float
            Scale-dependent growth factor D_cb(k,z)
        """
        self.compute(["fourier"])

        # if needed, convert k to units of 1/Mpc
        if h_units:
            k = k*self.ba.h

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

        # Get P_cb(k,z) and P_cb(k,0)
        if use_pk_lin == False:
            Pk_cb = self.pk_cb(k,z)
            Pk0_cb = self.pk_cb(k,0)
        else:
            Pk_cb = self.pk_cb_lin(k,z)
            Pk0_cb = self.pk_cb_lin(k,0)

        # Compute growth factor f
        D_cb = np.sqrt(Pk_cb/Pk0_cb)

        return D_cb

    def scale_independent_f_sigma8(self, z):
        """
        Return sigma8 * f(z)

        Return the scale-independent growth rate f(z) multiplied by
        sigma8(z) (dimensionless). The scale-independent growth rate
        f(z) is inferred from the scale-independent growth factor
        D(z), itself obtained by solving a background equation.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        scale_independent_f_sigma8 : float
            f(z)*sigma8(z) (single value or array)
        """
        return self.scale_independent_growth_factor_f(z)*self.sigma(8,z,h_units=True)

    def effective_f_sigma8(self, z, z_step=0.1):
        """
        Return sigma8(z) * f(k,z) estimated near k = 8 h/Mpc

        Return the growth rate f(k,z) estimated near k = 8 h/Mpc and multiplied by
        sigma8(z) (dimensionless). The product is actually directly inferred from sigma(R,z) as:

        sigma8(z) * f(k=8h/Mpc,z) = d sigma8 / d ln a = - (d sigma8 / dz)*(1+z)

        In the general case the quantity (d sigma8 / dz) is inferred
        from the two-sided finite difference between z+z_step and
        z-z_step. For z < z_step the step is reduced progressively
        down to z_step/10 while sticking to a double-sided
        derivative. For z< z_step/10 a single-sided derivative is used
        instead. (default: z_step=0.1)

        The input z can be a single value or an array of values. The
        output has the same shape.

        Parameters
        ----------
        z : float
            Redshift (single value or array)
        z_step : float, optional
            Step in redshift space for the finite difference method

        Returns
        -------
        effective_f_sigma8 : float
            sigma8(z) * f(k=8 h/Mpc,z) (single value or array)
        """
        self.compute(["fourier"])

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)
        if hasattr(z, "__len__"):
          out_array = np.empty_like(z,dtype=np.float64)
          for iz, redshift in enumerate(z):
            out_array[iz] = self.effective_f_sigma8(redshift, z_step=z_step)
          return out_array

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

    def effective_f_sigma8_spline(self, z, Nz=20):
        """
        Return sigma8(z) * f(k,z) estimated near k = 8 h/Mpc

        Return the growth rate f(k,z) estimated near k = 8 h/Mpc and multiplied by
        sigma8(z) (dimensionless). The product is actually directly inferred from sigma(R,z) as:

        sigma8(z) * f(k=8h/Mpc,z) = d sigma8 / d ln a = - (d sigma8 / dz)*(1+z)

        The derivative is inferred from a cubic spline method with Nz
        points, using the CubicSpline().derivative() method of the
        scipy.interpolate module (Nz=20 by default).

        The input z can be a single value or an array of values. The
        output has the same shape.

        Parameters
        ----------
        z : float
            Redshift (single value or array)
        Nz : float, optional
            Number of samples for the spline method

        Returns
        -------
        effective_f_sigma8 : float
            sigma8(z) * f(k=8 h/Mpc,z) (single value or array)
        """
        self.compute(["fourier"])

        if hasattr(z, "__len__"):
          out_array = np.empty_like(z,dtype=np.float64)
          for iz, redshift in enumerate(z):
            out_array[iz] = self.effective_f_sigma8_spline(redshift, Nz=Nz)
          return out_array

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)
        if self.fo.ln_tau_size>0:
          z_max = self.z_of_tau(np.exp(self.fo.ln_tau[0]))
        else:
          z_max = 0

        if (z<0) or (z>z_max):
            raise CosmoSevereError("You asked for effective_f_sigma8 at a redshift %e outside of the computed range [0,%e]"%(z,z_max))

        if (z<0.1):
            z_array = np.linspace(0, 0.2, num = Nz)
        elif (z<z_max-0.1):
            z_array = np.linspace(z-0.1, z+0.1, num = Nz)
        else:
            z_array = np.linspace(z_max-0.2, z_max, num = Nz)

        sig8_array = self.sigma(8,z_array,h_units=True)
        return -CubicSpline(z_array,sig8_array).derivative()(z)*(1+z)

    def z_of_tau(self, tau):
        """
        Return the redshift corresponding to a given conformal time.

        Return the redshift z (units of Mpc) corresponding to a given conformal time tau (units of Mpc).
        The user can pass a single tau or an array of tau. The shape of the output z is set accordingly.

        Parameters
        ----------
        tau : float
            Conformal time (single value or array)

        Returns
        -------
        z : float
            Redshift (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        tauarr = np.atleast_1d(tau).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        z = np.empty_like(tauarr)
        for itau, tauval in enumerate(tauarr):
          if background_at_tau(&self.ba,tauval,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          z[itau] = 1./pvecback[self.ba.index_bg_a]-1.

        free(pvecback)

        return (z[0] if np.isscalar(tau) else z)

    def Hubble(self, z):
        """
        Return Hubble(z)

        Return the Hubble rate H(z) - more precisely, the inverse Hubble radius H(z)/c (units 1/Mpc).
        The user can pass a single z or an array of z. The shape of the output z is set accordingly.

        Parameters
        ----------
        z : float
            Redshift single value or array)

        Returns
        -------
        H : float
            Hubble rate (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        H = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          H[iz] = pvecback[self.ba.index_bg_H]

        free(pvecback)

        return (H[0] if np.isscalar(z) else H)

    def Om_m(self, z):
        """
        Return the fractional density Omega of matter at z

        Return the fractional density Omega(z) of total matter at redshift z (dimensionless).
        The user can pass a single z or an array of z. The shape of the output is set accordingly.

        Parameters
        ----------
        z : float
            Redshift single value or array)

        Returns
        -------
        Om_m : float
            Omega_m(z) (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        Om_m = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          Om_m[iz] = pvecback[self.ba.index_bg_Omega_m]

        free(pvecback)

        return (Om_m[0] if np.isscalar(z) else Om_m)

    def Om_b(self, z):
        """
        Return the fractional density Omega of baryons at z

        Return the fractional density Omega(z) of baryons at redshift z (dimensionless).
        The user can pass a single z or an array of z. The shape of the output is set accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        Om_b : float
            Omega_b(z) (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        Om_b = np.empty_like(zarr)
        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          Om_b[iz] = pvecback[self.ba.index_bg_rho_b]/pvecback[self.ba.index_bg_rho_crit]

        free(pvecback)

        return (Om_b[0] if np.isscalar(z) else Om_b)

    def Om_cdm(self, z):
        """
        Return the fractional density Omega of CDM at z

        Return the fractional density Omega(z) of CDM at redshift z (dimensionless).
        The user can pass a single z or an array of z. The shape of the output is set accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        Om_cdm : float
            Omega_cdm(z) (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        Om_cdm = np.zeros_like(zarr)

        if self.ba.has_cdm == True:

          pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
          for iz, redshift in enumerate(zarr):

              if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                  free(pvecback) #manual free due to error
                  raise CosmoSevereError(self.ba.error_message)

              Om_cdm[iz] = pvecback[self.ba.index_bg_rho_cdm]/pvecback[self.ba.index_bg_rho_crit]

          free(pvecback)

        return (Om_cdm[0] if np.isscalar(z) else Om_cdm)

    def Om_ncdm(self, z):
        """
        Return the fractional density Omega of non-cold dark matter at z

        Return the fractional density Omega(z) of non-cold dark matter at redshift z (dimensionless).
        The user can pass a single z or an array of z. The shape of the output is set accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        Om_ncdm : float
            Omega_ncdm(z) (single value or array)
        """
        self.compute(["background"])

        cdef int last_index #junk
        cdef double * pvecback

        zarr = np.atleast_1d(z).astype(np.float64)

        Om_ncdm = np.zeros_like(zarr)

        if self.ba.has_ncdm == True:

            pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

            for iz, redshift in enumerate(zarr):
              if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
                  free(pvecback) #manual free due to error
                  raise CosmoSevereError(self.ba.error_message)

              rho_ncdm = 0.
              for n in range(self.ba.N_ncdm):
                  rho_ncdm += pvecback[self.ba.index_bg_rho_ncdm1+n]
              Om_ncdm[iz] = rho_ncdm/pvecback[self.ba.index_bg_rho_crit]

            free(pvecback)

        return (Om_ncdm[0] if np.isscalar(z) else Om_ncdm)

    def ionization_fraction(self, z):
        """
        Return the electron ionization fraction x_e(z)

        Return the electron ionization fraction x_e (dimensionless)
        for a given redshift z. CLASS sticks to the standard
        definition x_e = n_free_electrons / n_H, such that x_e can be
        bigger than one due to Helium.
        The user can pass a single z or an array of z. The shape of the output is set accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        xe : float
            Electron ionization fraction (single value or array)
        """
        self.compute(["thermodynamics"])

        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        zarr = np.atleast_1d(z).astype(np.float64)
        xe = np.empty_like(zarr)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              free(pvecthermo) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          if thermodynamics_at_z(&self.ba,&self.th,redshift,inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
              free(pvecback) #manual free due to error
              free(pvecthermo) #manual free due to error
              raise CosmoSevereError(self.th.error_message)

          xe[iz] = pvecthermo[self.th.index_th_xe]

        free(pvecback)
        free(pvecthermo)

        return (xe[0] if np.isscalar(z) else xe)

    def baryon_temperature(self, z):
        """
        Return the baryon temperature T_b(z)

        Return the baryon temperature T_b (units of K) for a given redshift z.
        The user can pass a single z or an array of z. The shape of the output is set accordingly.

        Parameters
        ----------
        z : float
            Redshift (single value or array)

        Returns
        -------
        Tb : float
            Baryon temperature (single value or array)
        """
        self.compute(["thermodynamics"])

        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        zarr = np.atleast_1d(z).astype(np.float64)
        Tb = np.empty_like(zarr)

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        for iz, redshift in enumerate(zarr):
          if background_at_z(&self.ba,redshift,long_info,inter_normal,&last_index,pvecback)==_FAILURE_:
              free(pvecback) #manual free due to error
              free(pvecthermo) #manual free due to error
              raise CosmoSevereError(self.ba.error_message)

          if thermodynamics_at_z(&self.ba,&self.th,redshift,inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
              free(pvecback) #manual free due to error
              free(pvecthermo) #manual free due to error
              raise CosmoSevereError(self.th.error_message)

          Tb[iz] = pvecthermo[self.th.index_th_Tb]

        free(pvecback)
        free(pvecthermo)

        return (Tb[0] if np.isscalar(z) else Tb)

    def get_background(self):
        """
        Return all background quantities

        Return a dictionary of background quantities at all times.
        The name and list of quantities in the returned dictionary are
        defined in CLASS, in background_output_titles() and
        background_output_data(). The keys of the dictionary refer to
        redshift 'z', proper time 'proper time [Gyr]', conformal time
        'conf. time [Mpc]', and many quantities such as the Hubble
        rate, distances, densities, pressures, or growth factors. For
        each key, the dictionary contains an array of values
        corresponding to each sampled value of time.

        This function works for whatever request in the 'output'
        field, and even if 'output' was not passed or left blank.

        Parameters
        ----------
        None

        Returns
        -------
        background : dict
            Dictionary of all background quantities at each time
        """
        self.compute(["background"])

        cdef char *titles
        cdef double* data
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if background_output_titles(&self.ba, titles)==_FAILURE_:
            free(titles) #manual free due to error
            raise CosmoSevereError(self.ba.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.ba.bt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if background_output_data(&self.ba, number_of_titles, data)==_FAILURE_:
            free(titles) #manual free due to error
            free(data) #manual free due to error
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
        Return all thermodynamics quantities

        Return a dictionary of thermodynbamics quantities at all
        times.  The name and list of quantities in the returned
        dictionary are defined in CLASS, in
        thermodynamics_output_titles() and
        thermodynamics_output_data(). The keys of the dictionary
        refer to the scale factor 'scale factor a', redshift 'z',
        conformal time 'conf. time [Mpc]', and many quantities such as
        the electron ionization fraction, scattering rates,
        temperatures, etc. For each key, the dictionary contains an
        array of values corresponding to each sampled value of time.

        This function works for whatever request in the 'output'
        field, and even if 'output' was not passed or left blank.

        Parameters
        ----------
        None

        Returns
        -------
        thermodynamics : dict
            Dictionary of all thermodynamics quantities at each time
        """
        self.compute(["thermodynamics"])

        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if thermodynamics_output_titles(&self.ba, &self.th, titles)==_FAILURE_:
            free(titles) #manual free due to error
            raise CosmoSevereError(self.th.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.th.tt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if thermodynamics_output_data(&self.ba, &self.th, number_of_titles, data)==_FAILURE_:
            free(titles) #manual free due to error
            free(data) #manual free due to error
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
        Return primordial spectra

        Return the primordial scalar spectrum, and possibly the
        primordial tensor spectrum if 'modes' includes 't'.  The name
        and list of quantities in the returned dictionary are defined
        in CLASS, in primordial_output_titles() and
        primordial_output_data().  The keys are 'k [1/Mpc]',
        'P_scalar(k)' and possibly 'P_tensor(k)'. For each key, the
        dictionary contains an array of values corresponding to each
        sampled wavenumber.

        This function requires that 'output' is set to something, e.g. 'tCl' or 'mPk'.

        Parameters
        ----------
        None

        Returns
        -------
        primordial : dict
            Dictionary of primordial spectra at each wavenumber
        """
        self.compute(["primordial"])

        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if primordial_output_titles(&self.pt, &self.pm, titles)==_FAILURE_:
            free(titles) #manual free due to error
            raise CosmoSevereError(self.pm.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pm.lnk_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if primordial_output_data(&self.pt, &self.pm, number_of_titles, data)==_FAILURE_:
            free(titles) #manual free due to error
            free(data) #manual free due to error
            raise CosmoSevereError(self.pm.error_message)

        primordial = {}

        for i in range(number_of_titles):
            primordial[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                primordial[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return primordial

    def get_perturbations(self, return_copy=True):
        """
        Return transfer function evolution for selected wavenumbers

        Return an array of dictionaries of transfer functions at all
        times (dimensionless). There is one dictionary for each
        requested mode in the list 'scalars', 'vectors', tensors', and
        for each wavenumbers requested by the user with the input
        parameter 'k_output_values' (units of 1/Mpc). For instance,
        perturbations['scalars'][0] could be one of the dictionaries.

        For each requested mode and wavenumber, the name and list of
        quantities in the returned dictionary are defined in CLASS, in
        perturbations_output_titles() and perturbations_output_data(),
        sticking to the 'class' format. The keys of the dictionary
        refer to wavnumbers 'k (h/Mpc)', density fluctuations 'd_g',
        etc., velocity perturbations 't_g', etc., and metric
        perturbations. For each key, the dictionary contains an array
        of values corresponding to each sampled value of time. An
        example of such array would be:

        perturbations['scalars'][0]['d_g'].

        Do not enable 'return_copy=False' unless you know exactly what
        you are doing.  This will mean that you get access to the
        direct C pointers inside CLASS.  That also means that if CLASS
        is deallocated, your perturbations array will become invalid.

        This function requires setting 'output' to somethings, for instance 'tCl' or 'mPk'.

        Parameters
        ----------
        return_copy : bool, optional
            Whether to return an exact copy of a memory area defined by C pointers inside CLASS

        Returns
        -------
        perturbations : dict
            Array of dictionaries of all transfer functions at each time
        """
        self.compute(["perturbations"])

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
                        tmpdict[names[i]] = (np.asarray(data_mv[:,i]).copy() if return_copy else np.asarray(data_mv[:,i]))
                    tmparray.append(tmpdict)
            perturbations[mode] = tmparray

        return perturbations

    def get_transfer(self, z=0., output_format='class'):
        """
        Return all scalar transfer functions at z

        Return a dictionary of transfer functions at all wavenumbers
        for scalar perturbations and at redshift z.  The name and list
        of quantities in the returned dictionary are defined in CLASS,
        in perturbations_output_titles() and
        perturbations_output_data(), using either the 'class' format
        (default) or 'camb' format. The keys of the dictionary refer
        to wavenumbers 'k (h/Mpc)', and in 'class' format, to density
        fluctuations 'd_g', etc., velocity perturbations 't_g', etc.,
        and metric perturbations. In 'camb' format, besides
        wavenumbers 'k (h/Mpc)', the returned transfer functions are
        '-T_g/k2', etc.  For each key, the dictionary contains an
        array of values corresponding to each sampled value of
        wavenumber.

        This function works if 'output' contains at least 'mTk' for
        density transfer function and/or 'vTk' for velocity transfer
        functions. To get the transfer functions at some z>0, the user
        must set 'z_pk' or 'z_max_pk' to a value equal or bigger to
        the requested one.

        Parameters
        ----------
        z : float, optional
            Redshift (0 by default)

        output_format : str, optional
            Output format, one of 'class' (defualt) or 'camb'

        Returns
        -------
        transfers : dict
            Dictionary of all transfer functions at each wavenumber
        """
        self.compute(["transfer"])

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
            free(titles) #manual free due to error
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
                free(titles) #manual free due to error
                free(data) #manual free due to error
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
        Return a dictionary of numbers computed internally by CLASS

        Return a dictionary containing an entry for each of the names defined in the input list 'names'. These names have to match the list defined inside this function, that is, belong to the list: 'h', 'H0' (units of km/s/Mpc), 'Omega0_lambda', 'Omega_Lambda', 'Omega0_fld', 'age' (units of Gyr), 'conformal_age' (units of Mpc), 'm_ncdm_in_eV' (units of eV), 'm_ncdm_tot' (units of eV), 'Neff', 'Omega_m', 'omega_m', 'xi_idr', 'N_dg', 'Gamma_0_nadm', 'a_dark', 'tau_reio', 'z_reio', 'z_rec', 'tau_rec' (units of Mpc), 'rs_rec' (units of Mpc), 'rs_rec_h' (units of Mpc/h), 'ds_rec' (units of Mpc), 'ds_rec_h' (units of Mpc/h), 'ra_rec' (units of Mpc), 'ra_rec_h' (units of Mpc/h), 'da_rec'(units of Mpc), 'da_rec_h' (units of Mpc/h), 'z_star', 'tau_star' (units of Mpc), 'rs_star' (units of Mpc), 'ds_star' (units of Mpc), 'ra_star' (units of Mpc), 'da_star' (units of Mpc), 'rd_star' (units of Mpc), 'z_d', 'tau_d' (units of Mpc), 'ds_d' (units of Mpc), 'ds_d_h' (units of Mpc/h), 'rs_d' (units of Mpc), 'rs_d_h' (units of Mpc/h), 'conf_time_reio' (units of Mpc), '100*theta_s', '100*theta_star', 'theta_s_100', 'theta_star_100', 'YHe', 'n_e', 'A_s', 'ln10^{10}A_s', 'ln_A_s_1e10', 'n_s', 'alpha_s', 'beta_s', 'r', 'r_0002', 'n_t', 'alpha_t', 'V_0', 'V_1', 'V_2', 'V_3', 'V_4', 'epsilon_V', 'eta_V', 'ksi_V^2', 'exp_m_2_tau_As', 'phi_min', 'phi_max', 'sigma8', 'sigma8_cb', 'k_eq' (units of 1/Mpc), 'a_eq', 'z_eq', 'H_eq' (units of 1/Mpc), 'tau_eq' (units of Mpc), 'g_sd', 'y_sd', 'mu_sd'.

        For instance, to get the age of the universe in Gyr and the Hubble parameter in km/s/Mpc, the user may call .get_current_derived_parameters(['age','H0']), store the result in a dictionary 'derived', and retrieve these values as age=derived['age'] and H0=derived['H0'].

        Parameters
        ----------
        names : str list
            Names defined inside this function, associated to quantities computed internally by CLASS

        Returns
        -------
        derived : dict
            Dictionary whose keys are the input names, and whose values are numbers extracted from the CLASS structures
        """
        if type(names) != type([]):
            raise TypeError("Deprecated")

        self.compute(["thermodynamics"])

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
            elif name == 'conf_time_reio':
                value = self.th.conf_time_reio
            elif name == '100*theta_s':
                value = 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)
            elif name == '100*theta_star':
                value = 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)
            elif name == 'theta_s_100':
                value = 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)
            elif name == 'theta_star_100':
                value = 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)
            elif name == 'YHe':
                value = self.th.YHe
            elif name == 'n_e':
                value = self.th.n_e
            elif name == 'A_s':
                value = self.pm.A_s
            elif name == 'ln10^{10}A_s':
                value = log(1.e10*self.pm.A_s)
            elif name == 'ln_A_s_1e10':
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
                self.compute(["fourier"])
                if (self.pt.has_pk_matter == _FALSE_):
                    raise CosmoSevereError("No power spectrum computed. In order to get sigma8, you must add mPk to the list of outputs.")
                value = self.fo.sigma8[self.fo.index_pk_m]
            elif name == 'sigma8_cb':
                self.compute(["fourier"])
                if (self.pt.has_pk_matter == _FALSE_):
                    raise CosmoSevereError("No power spectrum computed. In order to get sigma8_cb, you must add mPk to the list of outputs.")
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
                self.compute(["distortions"])
                if (self.sd.has_distortions == _FALSE_):
                    raise CosmoSevereError("No spectral distortions computed. In order to get g_sd, you must add sd to the list of outputs.")
                value = self.sd.sd_parameter_table[0]
            elif name == 'y_sd':
                self.compute(["distortions"])
                if (self.sd.has_distortions == _FALSE_):
                    raise CosmoSevereError("No spectral distortions computed. In order to get y_sd, you must add sd to the list of outputs.")
                value = self.sd.sd_parameter_table[1]
            elif name == 'mu_sd':
                self.compute(["distortions"])
                if (self.sd.has_distortions == _FALSE_):
                    raise CosmoSevereError("No spectral distortions computed. In order to get mu_sd, you must add sd to the list of outputs.")
                value = self.sd.sd_parameter_table[2]
            else:
                raise CosmoSevereError("%s was not recognized as a derived parameter" % name)
            derived[name] = value
        return derived

    def nonlinear_scale(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        Return the wavenumber associated to the nonlinearity scale for an array of redshifts

        Return the wavenumber k_nl(z) associated to the nonlinearity
        scale for an array of z_size redshifts. The nonlinearity
        scale is defined and computed within the Halofit or HMcode
        external modules.

        Parameters
        ----------
        z : float array
            Array of requested redshifts
        z_size : int
            Size of the redshift array

        Returns
        -------
        k_nl : float array
            Array of k_nl (z)
        """
        self.compute(["fourier"])

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
        Return the wavenumber associated to the nonlinearity scale (computed using the CDM+baryon spectrum) for an array of redshifts

        Return the wavenumber k_nl_cb(z) associated to the
        nonlinearity scale (computed using the CDM+baryon spectrum)
        for an array of z_size redshifts. The nonlinearity scale is
        defined and computed within the Halofit or HMcode external
        modules.

        Parameters
        ----------
        z : float array
            Array of requested redshifts
        z_size : int
            Size of the redshift array

        Returns
        -------
        k_nl_cb  : float array
            Array of k_nl_cb(z)
        """
        self.compute(["fourier"])

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
        """
        Return the total matter power spectrum P_m(k,z) computed with a fast method on a k and z array

        Return the total matter power spectrum P_m(k,z) (in Mpc**3)
        for an array of k values (in 1/Mpc) and z values, computed
        using a faster algorithm that with .get_pk(). The output is a
        flattened array index as pk[index_z*k_size+index_k].

        This function requires that the 'ouput' field contains at
        least 'mPk'.

        If 'non linear' is set to a non-linear method in CLASS, this
        function returns either the linear or non-linear power
        spectrum, depending on the argument 'nonlinear'. If 'non
        linear' is not set in CLASS, the argument 'nonlinear' should
        be set to 'False'.

        Parameters
        ----------
        k : float array
            Wavenumber array (one-dimensional, size k_size)
        z : float array
            Redshift array (one-dimensional, size z_size)
        k_size : int
            Number of k values
        z_size : int
            Number of redshift values
        nonlinear : bool
            Whether to return the nonlinear or linear power spectrum

        Returns
        -------
        pk : float array
            Flattened array of total matter power spectrum indexed as pk[index_z*k_size+index_k]
        """
        self.compute(["fourier"])
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == False:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk

    def get_pk_cb_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """
        Return the CDM+baryon power spectrum P_cb(k,z) computed with a fast method on a k and z array

        Return the CDM+baryon power spectrum P_cb(k,z) (in Mpc**3)
        for an array of k values (in 1/Mpc) and z values, computed
        using a faster algorithm that with .get_pk_cb(). The output is a
        flattened array index as pk_cb[index_z*k_size+index_k].

        This function requires that the 'ouput' field contains at
        least 'mPk'.

        If 'non linear' is set to a non-linear method in CLASS, this
        function returns either the linear or non-linear power
        spectrum, depending on the argument 'nonlinear'. If 'non
        linear' is not set in CLASS, the argument 'nonlinear' should
        be set to 'False'.

        Parameters
        ----------
        k : float array
            Wavenumber array (one-dimensional, size k_size)
        z : float array
            Redshift array (one-dimensional, size z_size)
        k_size : int
            Number of k values
        z_size : int
            Number of redshift values
        nonlinear : bool
            Whether to return the nonlinear or linear power spectrum

        Returns
        -------
        pk_cb : float array
            Flattened array of CDM+baryon power spectrum indexed as pk_cb[index_z*k_size+index_k]
        """
        self.compute(["fourier"])
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == False:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            fourier_pks_at_kvec_and_zvec(&self.ba, &self.fo, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk_cb

    def spectral_distortion_amplitudes(self):
        """
        Return the spectral distortion amplitudes g, mu, y, etc.

        Return the spectral distortion amplitude for each distorsion
        type g, mu, y, and residuals. The number of outputs beyond (g,
        mu, y) is the number of PCA components requested by the user
        using the input parameteer 'sd_PCA_size'

        This function requires that the 'ouput' field contains at least 'Sd'.

        Parameters
        ----------
        None

        Returns
        -------
        sd_type_amps : float array
            Spectral distortion amplitudes
        """
        self.compute(["distortions"])
        if self.sd.type_size == 0:
          raise CosmoSevereError("No spectral distortions have been calculated. Check that the output contains 'Sd' and the compute level is at least 'distortions'.")
        cdef np.ndarray[DTYPE_t, ndim=1] sd_type_amps = np.zeros(self.sd.type_size,'float64')
        for i in range(self.sd.type_size):
          sd_type_amps[i] = self.sd.sd_parameter_table[i]
        return sd_type_amps

    def spectral_distortion(self):
        """
        Return the shape of the total spectral distortion as a function of frequency

        Return the shape of the total spectral distortion (units of
        10^26 W m^-2 Hz^-1 sr^-1) as a function of frequency (units of
        GHz).

        This function requires that the 'ouput' field contains at least 'Sd'.

        Parameters
        ----------
        None

        Returns
        -------
        sd_nu : float array
            Array of frequencies nu/[1 GHz] (units of GHz)

        sd_amp : float array
            Array of total spectral distortion at these frequencies (units of 10^26 W m^-2 Hz^-1 sr^-1)
        """
        self.compute(["distortions"])
        if self.sd.x_size == 0:
          raise CosmoSevereError("No spectral distortions have been calculated. Check that the output contains 'Sd' and the compute level is at least 'distortions'.")
        cdef np.ndarray[DTYPE_t, ndim=1] sd_amp = np.zeros(self.sd.x_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sd_nu = np.zeros(self.sd.x_size,'float64')
        for i in range(self.sd.x_size):
          sd_amp[i] = self.sd.DI[i]*self.sd.DI_units*1.e26
          sd_nu[i] = self.sd.x[i]*self.sd.x_to_nu
        return sd_nu,sd_amp

    # Depracated functionalities
    def Omega0_m(self):
        """
        Return the Omega of total matter

        Return the fractional density Omega (dimensionless) of total
        non-relativistic matter (baryons, dark matter, massive
        neutrinos...) evaluated today. Strictly identical to the
        previously defined property Omega_m, but we leave it not to
        break compatibility.

        Parameters
        ----------
        None

        Returns
        -------
        Omega0_m : float
            Omega of total matter
        """
        return self.ba.Omega0_m


    def Omega0_k(self):
        """
        Return the Omega of curvature

        Return the effective fractional density Omega (dimensionless) of curvature evaluated today.
        Strictly identical to the previously defined
        property Omega_k, but we leave it not to break
        compatibility.

        Parameters
        ----------
        None

        Returns
        -------
        Omega0_k : float
            Omega of curvature
        """
        return self.ba.Omega0_k

    def Omega0_cdm(self):
        """
        Return the Omega of cdm

        Return the fractional density Omega (dimensionless) of CDM
        evaluated today. Strictly identical to the previously defined
        property Omega_cdm, but we leave it not to break
        compatibility.

        Parameters
        ----------
        None

        Returns
        -------
        Omega0_cdm : float
            Omega of CDM
        """
        return self.ba.Omega0_cdm

    # Source functions
    def get_sources(self):
        """
        Return the source functions stored in the perturbation module at all sampled (k, tau)

        After integrating all perturbations (more precisely, all
        transfer functions) over time, the CLASS perturbation module
        stores a list of the most important combinations of them on a
        grid of (k,tau) values. This includes the combinations used
        later to compute the CMB and LSS observables. This function
        retrieves these source functions directly from the CLASS
        structure 'perturbations', and returns them together with the
        values of k and tau at which they are sampled.

        The output is a dictionary whose keys belong to the list:
        -'p' for the CMB polarization source
        -'phi' for the metric fluctuation phi of the Newtonian gauge
        -'phi_plus_psi' for the metric fluctuation phi+psi of the Newtonian gauge
        -'phi_prime' for the metric fluctuation phi' of the Newtonian gauge
        -'psi' for the metric fluctuation psi of the Newtonian gauge
        -'H_T_Nb_prime' for the metric fluctuation H_T' of the Nboisson gauge
        -'k2gamma_Nb' for the k^2 gamma of the Nboisson gauge
        -'h' for the metric fluctuation h of the synchronouys gauge
        -'h_prime'  for the metric fluctuation h' of the synchronouys gauge
        -'eta' for the metric fluctuation eta of the synchronouys gauge
        -'eta_prime'  for the metric fluctuation eta' of the synchronouys gauge
        -'delta_tot'  for the total density fluctuation
        -'delta_m' for the non-relativistic density fluctuation
        -'delta_cb' for the CDM+baryon density fluctuation
        -'delta_g' for the photon density fluctuation
        -'delta_b' for the baryon density fluctuation
        -'delta_cdm' for the CDM density fluctuation
        -'delta_idm' for the interacting DM density fluctuation
        -'delta_dcdm' for the decaying DM density fluctuation
        -'delta_fld' for the fluid density fluctuation
        -'delta_scf' for the scalar field density fluctuation
        -'delta_dr' for the decay radiation density fluctuation
        -'delta_ur' for the ultra-relativistic density fluctuation
        -'delta_idr' for the interacting dark radiation density fluctuation
        -'delta_ncdm[i]' for the ncdm[i] density fluctuation
        -'theta_tot' for the total velocity divergence
        -'theta_m' for the non-relativistic velocity divergence
        -'theta_cb' for the CDM+baryon velocity divergence
        -'theta_g' for the photon velocity divergence
        -'theta_b' for the baryon velocity divergence
        -'theta_cdm' for the CDM velocity divergence
        -'theta_idm' for the interacting DM velocity divergence
        -'theta_dcdm' for the decaying DM velocity divergence
        -'theta_fld' for the fluid velocity divergence
        -'theta_scf' for the scalar field velocity divergence
        -'theta_dr' for the decayt radiation velocity divergence
        -'theta_ur' for the ultra-relativistic velocity divergence
        -'theta_ncdm[i]' for the ncdm[i] velocity divergence

        The list of actual keys in the returned dictionary depends on
        what the user requested in the 'output' field.

        For each key, the returned dictionary contains an array of sources, indexed as:

        sources['key'][index_k][index_tau]

        Parameters
        ----------
        None

        Returns
        -------
        sources : dict
            Dictionary containing the source functions sources['key'][index_k][index_tau] (dimensionless)
        numpy array : k_array
            Array of k values (units of 1/Mpc)
        numpy array : tau_array
            Array of tau values (units of Mpc)
        """
        self.compute(["fourier"])
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
            int tp_size = self.pt.tp_size[index_md];
            double *** sources_ptr = self.pt.sources;
            double [:,:] tmparray = np.zeros((k_size, tau_size)) ;
            double [:] k_array = np.zeros(k_size);
            double [:] tau_array = np.zeros(tau_size);

        names = []

        for index_k in range(k_size):
            k_array[index_k] = k[index_k]
        for index_tau in range(tau_size):
            tau_array[index_tau] = tau[index_tau]

        indices = []

        if self.pt.has_source_t:
            indices.extend([
                self.pt.index_tp_t0,
                self.pt.index_tp_t1,
                self.pt.index_tp_t2
                ])
            names.extend([
                "t0",
                "t1",
                "t2"
                ])
        if self.pt.has_source_p:
            indices.append(self.pt.index_tp_p)
            names.append("p")
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
        if self.pt.has_source_H_T_Nb_prime:
            indices.append(self.pt.index_tp_H_T_Nb_prime)
            names.append("H_T_Nb_prime")
        if self.pt.index_tp_k2gamma_Nb:
            indices.append(self.pt.index_tp_k2gamma_Nb)
            names.append("k2gamma_Nb")
        if self.pt.has_source_h:
            indices.append(self.pt.index_tp_h)
            names.append("h")
        if self.pt.has_source_h_prime:
            indices.append(self.pt.index_tp_h_prime)
            names.append("h_prime")
        if self.pt.has_source_eta:
            indices.append(self.pt.index_tp_eta)
            names.append("eta")
        if self.pt.has_source_eta_prime:
            indices.append(self.pt.index_tp_eta_prime)
            names.append("eta_prime")
        if self.pt.has_source_delta_tot:
            indices.append(self.pt.index_tp_delta_tot)
            names.append("delta_tot")
        if self.pt.has_source_delta_m:
            indices.append(self.pt.index_tp_delta_m)
            names.append("delta_m")
        if self.pt.has_source_delta_cb:
            indices.append(self.pt.index_tp_delta_cb)
            names.append("delta_cb")
        if self.pt.has_source_delta_g:
            indices.append(self.pt.index_tp_delta_g)
            names.append("delta_g")
        if self.pt.has_source_delta_b:
            indices.append(self.pt.index_tp_delta_b)
            names.append("delta_b")
        if self.pt.has_source_delta_cdm:
            indices.append(self.pt.index_tp_delta_cdm)
            names.append("delta_cdm")
        if self.pt.has_source_delta_idm:
            indices.append(self.pt.index_tp_delta_idm)
            names.append("delta_idm")
        if self.pt.has_source_delta_dcdm:
            indices.append(self.pt.index_tp_delta_dcdm)
            names.append("delta_dcdm")
        if self.pt.has_source_delta_fld:
            indices.append(self.pt.index_tp_delta_fld)
            names.append("delta_fld")
        if self.pt.has_source_delta_scf:
            indices.append(self.pt.index_tp_delta_scf)
            names.append("delta_scf")
        if self.pt.has_source_delta_dr:
            indices.append(self.pt.index_tp_delta_dr)
            names.append("delta_dr")
        if self.pt.has_source_delta_ur:
            indices.append(self.pt.index_tp_delta_ur)
            names.append("delta_ur")
        if self.pt.has_source_delta_idr:
            indices.append(self.pt.index_tp_delta_idr)
            names.append("delta_idr")
        if self.pt.has_source_delta_ncdm:
            for incdm in range(self.ba.N_ncdm):
              indices.append(self.pt.index_tp_delta_ncdm1+incdm)
              names.append("delta_ncdm[{}]".format(incdm))
        if self.pt.has_source_theta_tot:
            indices.append(self.pt.index_tp_theta_tot)
            names.append("theta_tot")
        if self.pt.has_source_theta_m:
            indices.append(self.pt.index_tp_theta_m)
            names.append("theta_m")
        if self.pt.has_source_theta_cb:
            indices.append(self.pt.index_tp_theta_cb)
            names.append("theta_cb")
        if self.pt.has_source_theta_g:
            indices.append(self.pt.index_tp_theta_g)
            names.append("theta_g")
        if self.pt.has_source_theta_b:
            indices.append(self.pt.index_tp_theta_b)
            names.append("theta_b")
        if self.pt.has_source_theta_cdm:
            indices.append(self.pt.index_tp_theta_cdm)
            names.append("theta_cdm")
        if self.pt.has_source_theta_idm:
            indices.append(self.pt.index_tp_theta_idm)
            names.append("theta_idm")
        if self.pt.has_source_theta_dcdm:
            indices.append(self.pt.index_tp_theta_dcdm)
            names.append("theta_dcdm")
        if self.pt.has_source_theta_fld:
            indices.append(self.pt.index_tp_theta_fld)
            names.append("theta_fld")
        if self.pt.has_source_theta_scf:
            indices.append(self.pt.index_tp_theta_scf)
            names.append("theta_scf")
        if self.pt.has_source_theta_dr:
            indices.append(self.pt.index_tp_theta_dr)
            names.append("theta_dr")
        if self.pt.has_source_theta_ur:
            indices.append(self.pt.index_tp_theta_ur)
            names.append("theta_ur")
        if self.pt.has_source_theta_idr:
            indices.append(self.pt.index_tp_theta_idr)
            names.append("theta_idr")
        if self.pt.has_source_theta_ncdm:
            for incdm in range(self.ba.N_ncdm):
              indices.append(self.pt.index_tp_theta_ncdm1+incdm)
              names.append("theta_ncdm[{}]".format(incdm))

        for index_type, name in zip(indices, names):
            tmparray = np.empty((k_size,tau_size))
            for index_k in range(k_size):
                for index_tau in range(tau_size):
                    tmparray[index_k][index_tau] = sources_ptr[index_md][index_ic*tp_size+index_type][index_tau*k_size + index_k];

            sources[name] = np.asarray(tmparray)

        return (sources, np.asarray(k_array), np.asarray(tau_array))
