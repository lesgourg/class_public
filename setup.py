from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import re

# get CLASS version from header
with open('include/common.h') as f:
    for line in f:
        match = re.search('#define _VERSION_ "v(.*)"', line)
        if match:
            CLASS_VERSION = match.group(1)
            break

# include directories for compilation
CLASSY_INCLUDE_DIRS = [
    numpy.get_include(),
    'include',
    'hyrec',
]

# source files for the classy extension
CLASSY_SOURCES = [
    'python/classy.pyx',
    'python/dir.c',
    'source/input.c',
    'source/background.c',
    'source/thermodynamics.c',
    'source/perturbations.c',
    'source/primordial.c',
    'source/nonlinear.c',
    'source/transfer.c',
    'source/spectra.c',
    'source/lensing.c',
    'tools/growTable.c',
    'tools/dei_rkck.c',
    'tools/sparse.c',
    'tools/evolver_rkck.c',
    'tools/evolver_ndf15.c',
    'tools/arrays.c',
    'tools/parser.c',
    'tools/quadrature.c',
    'tools/hyperspherical.c',
    'tools/common.c',
    'tools/trigonometric_integrals.c',
    'hyrec/hyrectools.c',
    'hyrec/helium.c',
    'hyrec/hydrogen.c',
    'hyrec/history.c',
]

# the setup script with one extension
setup(
    name='classy',
    version=CLASS_VERSION,
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    install_requires=[
        'numpy',
    ],
    ext_modules=cythonize([
        Extension(
            'classy',
            CLASSY_SOURCES,
            include_dirs=CLASSY_INCLUDE_DIRS,
            define_macros=[
                ('CLASSY_BUILD', None),
                ('__CLASSDIR__', 'get_classy_dir()'),
                ('HYREC', None),
            ],
        ),
    ]),
    packages=[
        'classy.bbn',
    ],
    package_dir={
        'classy': '',
    },
    package_data={
        'classy.bbn': ['*.dat'],
    },
)
