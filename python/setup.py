from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as nm
import os
import subprocess as sbp
import os.path as osp

gccpath = osp.normpath(osp.dirname(sbp.check_output("gcc -print-libgcc-file-name",shell=True)))
print gccpath

setup(
      name='classy',
      description='Python interface to the Cosmological Boltzmann code CLASS',
      url='http://www.class-code.net',
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("classy", ["classy.pyx"],
                               include_dirs = [nm.get_include(),"../include"],
                               libraries=["class"],library_dirs=["../",gccpath],
                   extra_link_args=['-lgomp'],
                               )],
      data_files=(('bbn',['../bbn/sBBN.dat']),)
)
