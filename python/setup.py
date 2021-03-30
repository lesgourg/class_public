from distutils.core import setup
#from distutils.extension import Extension
from Cython.Distutils import Extension
from Cython.Distutils import build_ext

import numpy as nm
import os
import subprocess as sbp
import os.path as osp

# Recover the gcc compiler
GCCPATH_STRING = sbp.Popen(
    ['gcc', '-print-libgcc-file-name'],
    stdout=sbp.PIPE).communicate()[0]
GCCPATH = osp.normpath(osp.dirname(GCCPATH_STRING)).decode()

liblist = ["class"]
MVEC_STRING = sbp.Popen(
    ['gcc', '-lmvec'],
    stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec","m"]

# define absolute paths
root_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
include_folder = os.path.join(root_folder, "include")
classy_folder = os.path.join(root_folder, "python")
heat_folder = os.path.join(os.path.join(root_folder, "external"),"heating")
recfast_folder = os.path.join(os.path.join(root_folder, "external"),"RecfastCLASS")
hyrec_folder = os.path.join(os.path.join(root_folder, "external"),"HyRec2020")

# Recover the CLASS version
with open(os.path.join(include_folder, 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            # get rid of the " and the v
            VERSION = line.split()[-1][2:-1]
            break

import sys

setup(
    name='classy',
    version=VERSION,
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    cmdclass={'build_ext': build_ext},
    package_dir={
        "classynet": os.path.join(classy_folder, "nn"),
        "classynet.generate": os.path.join(classy_folder, "nn", "generate"),
        "classynet.training": os.path.join(classy_folder, "nn", "training"),
        "classynet.models": os.path.join(classy_folder, "nn", "models"),
        "classynet.data_providers": os.path.join(classy_folder, "nn", "data_providers"),
        "classynet.testing": os.path.join(classy_folder, "nn", "testing"),
        "classynet.plotting": os.path.join(classy_folder, "nn", "plotting"),
        "classynet.tests": os.path.join(classy_folder, "nn", "tests"),
    },
    packages=[
        "classynet",
        "classynet.generate",
        "classynet.training",
        "classynet.models",
        "classynet.data_providers",
        "classynet.testing",
        "classynet.plotting",
        "classynet.tests",
    ],
    ext_modules=[
        Extension("classy", [os.path.join(classy_folder, "classy.pyx")],
                  include_dirs=[nm.get_include(), include_folder, heat_folder, recfast_folder, hyrec_folder],
                  libraries=liblist,
                  library_dirs=[root_folder, GCCPATH],
                  extra_link_args=['-lgomp'],
                  cython_directives={'language_level': "3" if sys.version_info.major>=3 else "2"}
                  ),
        Extension("classynet.lhs",
                  [os.path.join(classy_folder, "nn", "lhs_python.cpp")],
                  extra_compile_args=["-fopenmp"],
                  extra_link_args=["-fopenmp"],
                  include_dirs=[nm.get_include()],
                  )
    ],
    #data_files=[('bbn', ['../bbn/sBBN.dat'])]
)
