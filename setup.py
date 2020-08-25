from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build import build
from distutils.errors import CompileError
from Cython.Distutils import build_ext

import numpy as nm
import os
import subprocess

# build class that calls `make libclass.a`
class my_build(build):
    user_options = build.user_options + [
        ('ompflag=', None, 'compiler flag for OpenMP'),
    ]

    def initialize_options(self):
        build.initialize_options(self)
        self.ompflag = None

    def finalize_options(self):
        build.finalize_options(self)
        if "OMPFLAG" in os.environ:
            self.ompflag = os.environ.get('OMPFLAG')

    def run(self):
        make_libclass_a = ['make', 'libclass.a']
        if self.ompflag is not None:
            make_libclass_a += ['OMPFLAG='+self.ompflag]
        try:
            subprocess.check_call(make_libclass_a, stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise CompileError(e.output)
        build.run(self)

# Recover the CLASS version
with open(os.path.join('include', 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            # get rid of the " and the v
            VERSION = line.split()[-1][2:-1]
            break

# Define cython extension and fix Python version
classy_ext = Extension('classy', [os.path.join('python', 'classy.pyx')],
                           include_dirs=[nm.get_include(), 'include'],
                           libraries=['class'],
                           library_dirs=['.'])
import six
classy_ext.cython_directives = {'language_level': "3" if six.PY3 else "2"}
        
setup(
    name='classy',
    version=VERSION,
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    cmdclass={'build': my_build, 'build_ext': build_ext},
    ext_modules=[classy_ext],
    #data_files=[('bbn', ['../bbn/sBBN.dat'])]
)
