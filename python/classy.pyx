#########################################
# CLASS wrapper in python, 
# 13/02/2012
#
# Originally written by Karim Benabed
# Modified by Benjamin Audren, Julien Lesgourgues
#
# Note: be careful of the possible mix-up between Class and class. The first
# will always refer to the cosmological code, the second, always to the python
# object
#########################################

import numpy as np
import os
from math import exp,log
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
cimport cython 

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i
# Bunch of declarations from C to python. The idea here is to define only the
# quantities that will be used, for input, output or intermediate manipulation,
# by the python wrapper. For instance, in the precision structure, the only
# item used here is its error message. That is why nothing more is defined from
# this structure. The rest is internal in Class.
# If, for whatever reason, you need an other, existing parameter from Class,
# remember to add it inside this cdef.
cdef extern from "class.h":
  
  ctypedef char FileArg[40]
  
  ctypedef char* ErrorMsg
  
  cdef struct precision:
    ErrorMsg error_message

  cdef struct background  :
    ErrorMsg error_message 
    int bg_size
    int index_bg_ang_distance
    int index_bg_conf_distance
    int index_bg_H
    short long_info
    short inter_normal
    double T_cmb
    double h
    double age
    double conformal_age
    double * m_ncdm_in_eV
    double Neff
    double Omega0_b
    double Omega0_cdm
    double Omega0_ncdm_tot
    double Omega0_lambda
    double Omega0_fld    

  cdef struct thermo:
    ErrorMsg error_message 
    int th_size
    int index_th_xe
    int index_th_Tb
    short inter_normal
    double tau_reio
    double z_reio
    double z_rec
    double tau_rec
    double rs_rec
    double ds_rec
    double da_rec
    double YHe

  cdef struct perturbs      :
    ErrorMsg error_message 
    int has_pk_matter

  cdef struct bessels        :
    ErrorMsg error_message 

  cdef struct transfers          :
    ErrorMsg error_message 

  cdef struct primordial            :
    ErrorMsg error_message 
    double A_s
    double n_s
    double alpha_s
    double r
    double n_t
    double alpha_t   
    double V0
    double V1
    double V2
    double V3
    double V4

  cdef struct spectra              :
    ErrorMsg error_message 
    int l_max_tot
    int ln_k_size
    int ct_size
    int index_ct_tt
    int index_ct_te
    int index_ct_ee
    int index_ct_bb
    double* ln_k
    double sigma8

  cdef struct output                :
    ErrorMsg error_message 

  cdef struct lensing                  :
    int index_lt_tt
    int index_lt_te
    int index_lt_ee
    int index_lt_bb
    int has_lensed_cls
    int l_lensed_max
    ErrorMsg error_message 

  cdef struct nonlinear                    :
    int method
    ErrorMsg error_message 

  cdef struct file_content:
   char * filename
   int size
   FileArg * name
   FileArg * value
   short * read

  void lensing_free(        void*)
  void spectra_free(        void*)
  void primordial_free(     void*)
  void transfer_free(       void*)
  void perturb_free(        void*)
  void thermodynamics_free( void*)
  void background_free(     void*)
  void bessel_free(         void*)
  void nonlinear_free(void*)
  
  cdef int _FAILURE_
  cdef int _FALSE_
  cdef int _TRUE_
  
  int bessel_init(void*,void*)
  int input_init(void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 void*,
                 char*)
  int background_init(void*,void*)
  int thermodynamics_init(void*,void*,void*)
  int perturb_init(void*,void*,void*,void*)
  int transfer_init(void*,void*,void*,void*,void*,void*)
  int primordial_init(void*,void*,void*)
  int spectra_init(void*,void*,void*,void*,void*,void*)
  int nonlinear_init(void*,void*,void*,void*,void*,void*,void*,void*,void*)
  int lensing_init(void*,void*,void*,void*,void*)
  
  int background_tau_of_z(void* pba, double z,double* tau)
  int background_at_tau(void* pba, double tau, short return_format, short inter_mode, int * last_index, double *pvecback)
  int thermodynamics_at_z(void * pba, void * pth, double z, short inter_mode, int * last_index, double *pvecback, double *pvecthermo)
  int spectra_cl_at_l(void* psp,double l,double * cl,double * * cl_md,double * * cl_md_ic)
  int lensing_cl_at_l(void * ple,int l,double * cl_lensed)
  int spectra_pk_at_z(
      void * pba,
      void * psp,
      int mode, 
      double z,
      double * output_tot,
      double * output_ic
      )

  int spectra_pk_at_k_and_z(
    void* pba,
    void * ppm,
    void * psp,
    double k,
    double z,
    double * pk,
    double * pk_ic)

  int nonlinear_pk_at_z(
    void * pnl,
    double z,
    double * pz_density,
    double * pz_velocity,
    double * pz_cross,
    int * k_size_at_z)

  int nonlinear_pk_at_k_and_z(
    void * pnl,
    double k,
    double z,
    double * pz_density,
    double * pz_velocity,
    double * pz_cross,
    int * k_size_at_z)

  int nonlinear_k_nl_at_z(void* pnl, double z,double* k_nl)

  cdef enum linear_or_logarithmic :
            linear
            logarithmic
            

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, do the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class ClassError(Exception):
  def __init__(self,error_message,init=False):
    print error_message
    if init:
      raise AttributeError
    else:
      raise NameError

# The actual Class wrapping, the only class we will call from MontePython
# (indeed the only one we will import, with the command:
# from classy import Class
cdef class Class:
  # List of used structures
  cdef precision pr
  cdef   background ba
  cdef   thermo th
  cdef   perturbs pt 
  cdef   bessels bs
  cdef   transfers tr
  cdef   primordial pm
  cdef   spectra sp
  cdef   output op
  cdef   lensing le
  cdef   nonlinear nl
  cdef   file_content fc
  
  cdef object ready # Flag
  cdef object _pars # Dictionary of the parameters
  cdef object ncp   # Keeps track of the structures initialized, in view of cleaning.
  
  # Defining two new properties to recover, respectively, the parameters used
  # or the age (set after computation). Follow this syntax if you want to
  # access other quantities. Alternatively, you can also define a method, and
  # call it (see _T_cmb method, at the very bottom).
  property pars:
    def __get__(self):
      return self._pars
  property age:
    def __get__(self):
      return self._age()
  property Omega_m:
    def __get__(self):
      return self._Omega_m()
  property Omega_nu:
    def __get__(self):
      return self.ba.Omega0_ncdm_tot
  property nonlinear_method:
    def __get__(self):
      return self.nl.method  

  def __init__(self,default=False):
    cdef char* dumc
    self.ready = False
    self._pars = {}
    self.fc.size=0
    self.fc.filename = <char*>malloc(sizeof(char)*30)
    assert(self.fc.filename!=NULL)
    dumc = "NOFILE"
    sprintf(self.fc.filename,"%s",dumc)
    self.ncp = set()
    
  # Set up the dictionary
  def set(self,*pars,**kars):
    if len(pars)==1:
      self._pars.update(dict(pars[0]))
    elif len(pars)!=0:
      raise ClassError("bad call")
    self._pars.update(kars)
    self.ready=False
  
  def empty(self):
    self._pars = {}
    self.ready=False
    
  def cleanup(self):
    if self.ready==False:
      return 
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
  def _struct_cleanup(self,ncp):
    if "lensing" in ncp:
      lensing_free(&self.le)
    if "nonlinear" in ncp:
      nonlinear_free (&self.nl)
    if "spectra" in ncp:
      spectra_free(&self.sp)
    if "primordial" in ncp:
      primordial_free(&self.pm)
    if "transfer" in ncp:
      transfer_free(&self.tr)
    if "perturb" in ncp:
      perturb_free(&self.pt)
    if "thermodynamics" in ncp:
      thermodynamics_free(&self.th)
    if "background" in ncp:
      background_free(&self.ba)
    if "bessel" in ncp:
      bessel_free(&self.bs)

  # Ensure the full module dependency
  def _check_task_dependency(self,ilvl):
    lvl = ilvl.copy()
    #print "before",lvl
    if "lensing" in lvl:
      lvl.add("nonlinear")
    if "nonlinear" in lvl:
      lvl.add("spectra")
    if "spectra" in lvl:
      lvl.add("primordial")
    if "primordial" in lvl:
      lvl.add("transfer")
    if "transfer" in lvl:
      lvl.add("bessel")
      lvl.add("perturb")
    if "perturb" in lvl:
      lvl.add("thermodynamics")
    if "thermodynamics" in lvl:
      lvl.add("background")
    if len(lvl)!=0 :
      lvl.add("input")
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
  def _compute(self,lvl=("lensing")):
    cdef ErrorMsg errmsg
    cdef int ierr
    cdef char* dumc
    
    lvl = self._check_task_dependency(set(lvl))
    
    if self.ready and self.ncp.issuperset(lvl):
      return
    self.ready = True
    
    self._fillparfile()
      
    # empty all
    #self._struct_cleanup(self.ncp)
    self.ncp=set()
    
    # compute
    if "input" in lvl:
      ierr = input_init(&self.fc,
                          &self.pr,
                          &self.ba,
                          &self.th,
                          &self.pt,
                          &self.bs,
                          &self.tr,
                          &self.pm,
                          &self.sp,
                          &self.nl,
                          &self.le,
                          &self.op,
                          errmsg)
      if ierr==_FAILURE_:
        raise ClassError(errmsg)
      self.ncp.add("input")      

      for i in range(self.fc.size):
        if self.fc.read[i] == _FALSE_:
          raise ClassError("Class did not read input parameter %s\n" % self.fc.name[i],init=True)
    
    if "background" in lvl:
      if background_init(&(self.pr),&(self.ba)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.ba.error_message)
        raise ClassError(self.ba.error_message)
      self.ncp.add("background") 
    
    if "thermodynamics" in lvl:
      if thermodynamics_init(&(self.pr),&(self.ba),&(self.th)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.th.error_message)
        raise ClassError(self.th.error_message)
      self.ncp.add("thermodynamics") 
  
    if "perturb" in lvl:
      if perturb_init(&(self.pr),&(self.ba),&(self.th),&(self.pt)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.pt.error_message)
        raise ClassError(self.pt.error_message)
      self.ncp.add("perturb") 
      
    if "bessel" in lvl:
      if bessel_init(&(self.pr),&(self.bs)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.bs.error_message)
        raise ClassError(self.bs.error_message)
      self.ncp.add("bessel") 
      
    if "transfer" in lvl:
      if transfer_init(&(self.pr),&(self.ba),&(self.th),&(self.pt),&(self.bs),&(self.tr)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.tr.error_message)
        raise ClassError(self.tr.error_message)
      self.ncp.add("transfer") 
      
    if "primordial" in lvl:
      if primordial_init(&(self.pr),&(self.pt),&(self.pm)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.pm.error_message)
        raise ClassError(self.pm.error_message)
      self.ncp.add("primordial") 
      
    if "spectra" in lvl:
      if spectra_init(&(self.pr),&(self.ba),&(self.pt),&(self.tr),&(self.pm),&(self.sp)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.sp.error_message)
        raise ClassError(self.sp.error_message)
      self.ncp.add("spectra")       

    if "nonlinear" in lvl:
      if (nonlinear_init(&self.pr,&self.ba,&self.th,&self.pt,&self.bs,&self.tr,&self.pm,&self.sp,&self.nl) == _FAILURE_):
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.nl.error_message)
        raise ClassError(self.nl.error_message)
      self.ncp.add("nonlinear") 
       
    if "lensing" in lvl:
      if lensing_init(&(self.pr),&(self.pt),&(self.sp),&(self.nl),&(self.le)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        #fprintf(stderr,"%s\n",self.le.error_message)
        raise ClassError(self.le.error_message)
      self.ncp.add("lensing") 
      
    # At this point, the cosmological instance contains everything needed. The
    # following functions are only to output the desired numbers
    return

  def raw_cl(self, lmax=-1,nofail=False):
    cdef int lmaxR 
    cdef double *lcl = <double*> calloc(self.sp.ct_size,sizeof(double))
    
    lmaxR = self.sp.l_max_tot
    if lmax==-1:
      lmax=lmaxR
    if lmax>lmaxR:
      if nofail:
        self._pars_check("l_max_scalar",lmax)
        self._compute(["lensing"])
      else:
        raise ClassError("Can only compute up to lmax=%d"%lmaxR)

    cl = {}
    for elem in ['tt','te','ee','bb']:
      cl[elem] = np.ndarray(lmax+1, dtype=np.double)
      cl[elem][:2]=0
    for ell from 2<=ell<lmax+1:
      if spectra_cl_at_l(&self.sp,ell,lcl,NULL,NULL) == _FAILURE_:
        raise ClassError(self.sp.error_message) 
      cl['tt'][ell] = lcl[self.sp.index_ct_tt]
      cl['te'][ell] = lcl[self.sp.index_ct_te]
      cl['ee'][ell] = lcl[self.sp.index_ct_ee]
      cl['bb'][ell] = lcl[self.sp.index_ct_bb]

    free(lcl)  
    return cl

  def lensed_cl(self, lmax=-1,nofail=False):
    cdef int lmaxR 
    cdef double *lcl = <double*> calloc(self.sp.ct_size,sizeof(double))
    lmaxR = self.le.l_lensed_max
    
    if lmax==-1:
      lmax=lmaxR
    if lmax>lmaxR:
      if nofail:
        self._pars_check("l_max_scalar",lmax)
        self._compute(["lensing"])
      else:
        raise ClassError("Can only compute up to lmax=%d"%lmaxR)
    
    cl = {}
    for elem in ['tt','te','ee','bb']:
      cl[elem] = np.ndarray(lmax+1, dtype=np.double)
      cl[elem][:2]=0
    for ell from 2<=ell<lmax+1:
      if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
        raise ClassError(self.le.error_message) 
      cl['tt'][ell] = lcl[self.le.index_lt_tt]
      cl['te'][ell] = lcl[self.le.index_lt_te]
      cl['ee'][ell] = lcl[self.le.index_lt_ee]
      cl['bb'][ell] = lcl[self.le.index_lt_bb]

    free(lcl)  
    return cl
    
  def z_of_r (self,z_array):
    cdef double tau=0.0
    cdef int last_index=0 #junk
    cdef double * pvecback
    r    = np.zeros(len(z_array),'float64')
    dzdr = np.zeros(len(z_array),'float64')

    pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

    i = 0
    for redshift in z_array:
      if background_tau_of_z(&self.ba,redshift,&tau)==_FAILURE_:
        raise ClassError(self.ba.error_message)

      if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
        raise ClassError(self.ba.error_message)

      # store r
      r[i] = pvecback[self.ba.index_bg_conf_distance]
      # store dz/dr = H
      dzdr[i] = pvecback[self.ba.index_bg_H]

      i += 1

    free(pvecback)
    return r[:],dzdr[:]

  # Gives the pk for a given (k,z)
  def _pk(self,double k,double z):
    cdef double pk
    cdef double pk_velo
    cdef double pk_cross
    cdef int dummy

    if (self.nl.method == 0):
       if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,k,z,&pk,NULL)==_FAILURE_:
         raise ClassError(self.sp.error_message)
    else:
       if nonlinear_pk_at_k_and_z(&self.nl,k,z,&pk,&pk_velo,&pk_cross,&dummy)==_FAILURE_:
          raise ClassError(self.nl.error_message)
    #free(junk)
    return pk

  def _get_pk(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
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
    return index
        
  def _age(self):
    self._compute(["background"])
    return self.ba.age
    
  def _h(self):
    return self.ba.h

  def _n_s(self):
    return self.pm.n_s  

  def _Omega_m(self):
    return self.ba.Omega0_b+self.ba.Omega0_cdm

  def _sigma8(self):
    self._compute(["spectra"])
    return self.sp.sigma8   

  def _angular_distance(self, z):
    cdef double tau
    cdef int last_index #junk
    cdef double * pvecback
    #D_A = np.zeros(np.shape(z_array)[0],'float64')

    pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

    if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
      raise ClassError(self.ba.error_message)

    if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
      raise ClassError(self.ba.error_message)

    D_A = pvecback[self.ba.index_bg_ang_distance]

    free(pvecback)
      
    return D_A

  def _ionization_fraction(self, z):
    cdef double tau
    cdef int last_index #junk
    cdef double * pvecback
    cdef double * pvecthermo

    pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
    pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

    if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
      raise ClassError(self.ba.error_message)

    if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
      raise ClassError(self.ba.error_message)

    if thermodynamics_at_z(&self.ba,&self.th,z,self.th.inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
      raise ClassError(self.th.error_message)

    xe = pvecthermo[self.th.index_th_xe] 

    free(pvecback)
    free(pvecthermo)

    return xe

  def _baryon_temperature(self, z):
    cdef double tau
    cdef int last_index #junk
    cdef double * pvecback
    cdef double * pvecthermo

    pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
    pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

    if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
      raise ClassError(self.ba.error_message)

    if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
      raise ClassError(self.ba.error_message)

    if thermodynamics_at_z(&self.ba,&self.th,z,self.th.inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
      raise ClassError(self.th.error_message)

    Tb = pvecthermo[self.th.index_th_Tb] 

    free(pvecback)
    free(pvecthermo)

    return Tb

  def _T_cmb(self):
    return self.ba.T_cmb

  def _Omega0_m(self):
    return self.ba.Omega0_b+self.ba.Omega0_cdm
  
  def get_current_derived_parameters(self,data):
    for elem in data.get_mcmc_parameters(['derived']):
      if elem == 'h':
        data.mcmc_parameters[elem]['current'] = self.ba.h
      elif elem == 'Omega0_lambda' or elem == 'Omega_Lambda':
        data.mcmc_parameters[elem]['current'] = self.ba.Omega0_lambda
      elif elem == 'Omega0_fld':
        data.mcmc_parameters[elem]['current'] = self.ba.Omega0_fld
      elif elem == 'age':
        data.mcmc_parameters[elem]['current'] = self.ba.age
      elif elem == 'conformal_age':
        data.mcmc_parameters[elem]['current'] = self.ba.conformal_age
      elif elem == 'm_ncdm_in_eV':
        data.mcmc_parameters[elem]['current'] = self.ba.m_ncdm_in_eV[0]
      elif elem == 'm_ncdm_tot':
        data.mcmc_parameters[elem]['current'] = self.ba.Omega0_ncdm_tot*self.ba.h*self.ba.h*93.14
      elif elem == 'Neff':
        data.mcmc_parameters[elem]['current'] = self.ba.Neff
      elif elem == 'tau_reio':
        data.mcmc_parameters[elem]['current'] = self.th.tau_reio
      elif elem == 'z_reio':
        data.mcmc_parameters[elem]['current'] = self.th.z_reio
      elif elem == 'z_rec':
        data.mcmc_parameters[elem]['current'] = self.th.z_rec
      elif elem == 'tau_rec':
        data.mcmc_parameters[elem]['current'] = self.th.tau_rec
      elif elem == 'rs_rec':
        data.mcmc_parameters[elem]['current'] = self.th.rs_rec
      elif elem == 'rs_rec_h':
        data.mcmc_parameters[elem]['current'] = self.th.rs_rec*self.ba.h
      elif elem == 'ds_rec':
        data.mcmc_parameters[elem]['current'] = self.th.ds_rec
      elif elem == 'ds_rec_h':
        data.mcmc_parameters[elem]['current'] = self.th.ds_rec*self.ba.h
      elif elem == 'ra_rec':
        data.mcmc_parameters[elem]['current'] = self.th.da_rec*(1.+self.th.z_rec)
      elif elem == 'ra_rec_h':
        data.mcmc_parameters[elem]['current'] = self.th.da_rec*(1.+self.th.z_rec)*self.ba.h
      elif elem == 'da_rec':
        data.mcmc_parameters[elem]['current'] = self.th.da_rec
      elif elem == 'da_rec_h':
        data.mcmc_parameters[elem]['current'] = self.th.da_rec*self.ba.h
      elif elem == 'YHe':
        data.mcmc_parameters[elem]['current'] = self.th.YHe
      elif elem == 'ne':
        data.mcmc_parameters[elem]['current'] = self.th.ne
      elif elem == 'A_s':
        data.mcmc_parameters[elem]['current'] = self.pm.A_s
      elif elem == 'ln10^{10}A_s':
        data.mcmc_parameters[elem]['current'] = log(1.e10*self.pm.A_s)
      elif elem == 'n_s':
        data.mcmc_parameters[elem]['current'] = self.pm.n_s
      elif elem == 'alpha_s':
        data.mcmc_parameters[elem]['current'] = self.pm.alpha_s
      elif elem == 'r':
        data.mcmc_parameters[elem]['current'] = self.pm.r
      elif elem == 'n_t':
        data.mcmc_parameters[elem]['current'] = self.pm.n_t
      elif elem == 'alpha_t':
        data.mcmc_parameters[elem]['current'] = self.pm.alpha_t
      elif elem == 'V_0':
        data.mcmc_parameters[elem]['current'] = self.pm.V0
      elif elem == 'V_1':
        data.mcmc_parameters[elem]['current'] = self.pm.V1
      elif elem == 'V_2':
        data.mcmc_parameters[elem]['current'] = self.pm.V2
      elif elem == 'V_3':
        data.mcmc_parameters[elem]['current'] = self.pm.V3
      elif elem == 'V_4':
        data.mcmc_parameters[elem]['current'] = self.pm.V4
      elif elem == 'exp_m_2_tau_As':
        data.mcmc_parameters[elem]['current'] = exp(-2.*self.th.tau_reio)*self.pm.A_s
      elif elem == 'sigma8':
        data.mcmc_parameters[elem]['current'] = self.sp.sigma8
      else:
        raise ClassError("%s was not recognized as a derived parameter" % elem)
    return

  def nonlinear_scale(self,np.ndarray[DTYPE_t,ndim=1] z,int z_size):
    cdef int index_z
    cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
    #cdef double *k_nl

    #k_nl = <double*> calloc(z_size,sizeof(double))
    for index_z in range(z_size):
      if nonlinear_k_nl_at_z(&self.nl,z[index_z],&k_nl[index_z]) == _FAILURE_:
        raise ClassError(self.nl.error_message)

    return k_nl
