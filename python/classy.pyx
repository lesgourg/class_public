import numpy as nm
cimport numpy as nm
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
cimport cython 

cdef extern from "class.h":
  
  ctypedef char FileArg[40]
  
  ctypedef char* ErrorMsg
  
  cdef struct precision:
    ErrorMsg error_message

  cdef struct background  :
    ErrorMsg error_message 
    double age
    
  cdef struct thermo:
    ErrorMsg error_message 

  cdef struct perturbs      :
    ErrorMsg error_message 
    int has_pk_matter

  cdef struct bessels        :
    ErrorMsg error_message 

  cdef struct transfers          :
    ErrorMsg error_message 

  cdef struct primordial            :
    ErrorMsg error_message 

  cdef struct spectra              :
    ErrorMsg error_message 
    int l_max_tot
    int ln_k_size
    double* ln_k

  cdef struct output                :
    ErrorMsg error_message 

  cdef struct lensing                  :
    int has_lensed_cls
    int l_lensed_max
    ErrorMsg error_message 

  cdef struct nonlinear                    :
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
  int nonlinear_init(void*,void*,void*,void*,void*,void*)
  int lensing_init(void*,void*,void*,void*,void*)
  
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
  cdef enum linear_or_logarithmic :
            linear
            logarithmic
            
class ClassError(Exception):
  pass

cdef class Class:
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
  
  cdef object ready
  cdef object _pars
  cdef object ncp
  
  property pars:
    def __get__(self):
      return self._pars
  property age:
    def __get__(self):
      return self._age()
      
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
    if default: self.set_default()
    
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
  
  def set_default(self):
    self._pars = {
                  "output":"tCl mPk",
                  "background_verbose" : 1,
                  "thermodynamics_verbose" : 1,
                  "perturbations_verbose" : 1,
                  "bessels_verbose" : 1,
                  "transfer_verbose" : 1,
                  "primordial_verbose" : 1,
                  "spectra_verbose" : 1,
                  "nonlinear_verbose" : 1,
                  "lensing_verbose" : 1,
                  "output_verbose": 1,
                  }
    self.ready=False

  
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
      i+=1
      
      
  def _struct_cleanup(self,ncp):
    #print "clean",ncp
    if "lensing" in ncp:
      #print "lensing"
      lensing_free(&self.le)
    if "nonlinear" in ncp:
      #print "nonlinear"
      nonlinear_free (&self.nl)
    if "spectra" in ncp:
      #print "spectra"
      spectra_free(&self.sp)
    if "primordial" in ncp:
      #print "primordial"
      primordial_free(&self.pm)
    if "transfer" in ncp:
      #print "transfer"
      transfer_free(&self.tr)
    if "perturb" in ncp:
      #print "perturb"
      perturb_free(&self.pt)
    if "thermodynamics" in ncp:
      #print "thermodynamics"
      thermodynamics_free(&self.th)
    if "background" in ncp:
      #print "background"
      background_free(&self.ba)
    if "bessel" in ncp:
      #print "bessel"
      bessel_free(&self.bs)

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
    #print "after",lvl
    return lvl
    
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
    
    if "background" in lvl:
      if background_init(&(self.pr),&(self.ba)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.ba.error_message)
      self.ncp.add("background") 
    
    if "thermodynamics" in lvl:
      if thermodynamics_init(&(self.pr),&(self.ba),&(self.th)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("thermodynamics") 
  
    if "perturb" in lvl:
      if perturb_init(&(self.pr),&(self.ba),&(self.th),&(self.pt)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("perturb") 
      
    if "bessel" in lvl:
      if bessel_init(&(self.pr),&(self.bs)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("bessel") 
      
    if "transfer" in lvl:
      if transfer_init(&(self.pr),&(self.ba),&(self.th),&(self.pt),&(self.bs),&(self.tr)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("transfer") 
      
    if "primordial" in lvl:
      if primordial_init(&(self.pr),&(self.pt),&(self.pm)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("primordial") 
      
    if "spectra" in lvl:
      if spectra_init(&(self.pr),&(self.ba),&(self.pt),&(self.tr),&(self.pm),&(self.sp)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("spectra")       

    if "nonlinear" in lvl:
      if (nonlinear_init(&self.pr,&self.ba,&self.th,&self.pm,&self.sp,&self.nl) == _FAILURE_):
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("nonlinear") 
       
    if "lensing" in lvl:
      if lensing_init(&(self.pr),&(self.pt),&(self.sp),&(self.nl),&(self.le)) == _FAILURE_:
        self._struct_cleanup(self.ncp)
        raise ClassError(self.th.error_message)
      self.ncp.add("lensing") 
      
    return

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
    
  def raw_cl(self, lmax=-1,nofail=False):
    cdef int lmaxR 
    cdef nm.ndarray cl
    cdef double lcl[4]
    
    if nofail:
      self._pars_check("output","tCl",True,True)
    
    self._compute(["input"])
    lmaxR = self.sp.l_max_tot
    if lmax==-1:
      lmax=lmaxR
    if lmax>lmaxR:
      if nofail:
        self._pars_check("l_max_scalar",lmax)
      else:
        raise ClassError("Can only compute up to lmax=%d"%lmaxR)
    self._compute(["spectra"])
    
    cl = nm.ndarray([4,lmax+1], dtype=nm.double)
    cl[:2]=0
    lcl[0]=lcl[1]=lcl[2]=lcl[3] = 0
    for ell from 2<=ell<lmax+1:
      if spectra_cl_at_l(&self.sp,ell,lcl,NULL,NULL) == _FAILURE_:
        raise ClassError(self.sp.error_message) 
      for md from 0<=md<4:
        cl[md,ell] = lcl[md]
    return cl

  def lensed_cl(self, lmax=-1,nofail=False):
    cdef int lmaxR 
    cdef nm.ndarray cl
    cdef double lcl[8]
    
    if nofail:
      self._pars_check("output","tCl",True,True)
      self._pars_check("output","lCl",True,True)
      self._pars_check("lensing","yes",False,True)

    self._compute(["input"])
    if self.le.has_lensed_cls==0:
      raise ClassError("No lensing effect computed")
    self._compute(["lensing"])
    lmaxR = self.le.l_lensed_max
    
    if lmax==-1:
      lmax=lmaxR
    if lmax>lmaxR:
      if nofail:
        self._pars_check("l_max_scalar",lmax)
        self._compute(["lensing"])
      else:
        raise ClassError("Can only compute up to lmax=%d"%lmaxR)
    
    cl = nm.ndarray([4,lmax+1], dtype=nm.double)
    cl[:2]=0
    lcl[0]=lcl[1]=lcl[2]=lcl[3] = 0
    for ell from 2<=ell<lmax+1:
      if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
        raise ClassError(self.le.error_message) 
      for md from 0<=md<4:
        cl[md,ell] = lcl[md]
      
    self._struct_cleanup(self.ncp)
    return cl
    
  def pk_l (self,double z=0,k=None,nofail=False):
    cdef nm.ndarray _k
    cdef nm.ndarray pk
    cdef double mpk,mk
    
    if nofail:
      self._pars_check("output","mPk",True,True)
      self._pars_check("z_pk",str(z),True,",")
      
    self._compute(["input"])
    if self.pt.has_pk_matter==0:
      raise ClassError("no pK computed")
    self._compute(["spectra"])

    if k==None:
      _k = nm.ndarray([self.sp.ln_k_size], dtype=nm.double)
      memcpy(nm.PyArray_DATA(_k),self.sp.ln_k,sizeof(double)*self.sp.ln_k_size)
      k = nm.exp(_k)
      pk = nm.ndarray([self.sp.ln_k_size], dtype=nm.double)
    
      if spectra_pk_at_z(&self.ba,&self.sp,linear,z,<double*>nm.PyArray_DATA(pk),NULL)==_FAILURE_:
        raise ClassError(self.sp.error_message)
      return nm.array((k,pk))
    
    else:
      pk = nm.ndarray([len(k)], dtype=nm.double)
      
      for i from 0<=i<len(k):
        mk = k[i]
        if spectra_pk_at_k_and_z(&self.ba,&self.pm,&self.sp,mk,z,&mpk,NULL)==_FAILURE_:
          raise ClassError(self.sp.error_message)
        pk[i] = mpk
      
      return nm.array((k,pk))
        
  def _age(self):
    self._compute(["background"])
    return self.ba.age
    
