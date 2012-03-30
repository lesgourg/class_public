from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as nm
import os

classlib = ["background.c","bessel.c","input.c","lensing.c","nonlinear.c","output.c","perturbations.c","primordial.c","spectra.c","thermodynamics.c","transfer.c","trg.c"]
classtool = ["arrays.c","dei_rkck.c","evolver_rkck.c","growTable.c","parser.c","quadrature.c","sparse.c"]

classlib = ["../source/"+src for src in classlib]
classtool = ["../tools/"+src for src in classtool]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("classy", ["classy.pyx"],
                             include_dirs = [nm.get_include(),"../include"],
                             libraries=["class"],library_dirs=["../","/opt/local/lib/gcc44/"],
			     extra_link_args=['-lgomp'],
                             )],
)
