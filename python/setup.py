from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy as nm
import os
import subprocess as sbp
import os.path as osp

# Recover the gcc compiler
GCCPATH_STRING = sbp.Popen(
    ['gcc', '-print-libgcc-file-name'],
    stdout=sbp.PIPE).communicate()[0]
GCCPATH = osp.normpath(osp.dirname(GCCPATH_STRING)).decode()

compile_args = ['-g', '-std=c++11']

liblist = ["class"]
MVEC_STRING = sbp.Popen(
    ['gcc', '-lmvec'],
    stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec","m"]

# define absolute paths
root_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
include_folders = [nm.get_include()]
include_folders.append(os.path.join(root_folder, "include"))
include_folders.append(os.path.join(root_folder, "source"))
include_folders.append(os.path.join(root_folder, "tools"))
include_folders.append(os.path.join(root_folder, "main"))
classy_folder = os.path.join(root_folder, "python")

# Recover the CLASS version
with open(os.path.join(include_folders[1], 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            # get rid of the " and the v
            VERSION = line.split()[-1][2:-1]
            break

# Define cython extension and fix Python version
classy_ext = Extension("classy", [os.path.join(classy_folder, "classy.pyx")],
                           include_dirs=include_folders,
                           libraries=liblist,
                           library_dirs=[root_folder, GCCPATH],
                           language="c++",
                           extra_compile_args=compile_args,
                           define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],)
import six

setup(
    name='classy',
    version=VERSION,
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    ext_modules=cythonize(
        classy_ext,
        language_level=3 if six.PY3 else 2,
        annotate=False,
    ),
)
