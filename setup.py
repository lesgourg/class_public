from setuptools import setup, Extension
from Cython.Distutils import build_ext

import numpy as np
import os
import subprocess as sbp
import sys
from pip._internal.commands.install import decide_user_install

def binaries_directory():
    import site
    """Return the installation directory, or None"""
    user_install = decide_user_install(use_user_site = True if ('--user' in sys.argv) else None, prefix_path = any(['--prefix' in x for x in sys.argv]), target_dir=any(['--target' in x for x in sys.argv]), root_path = None, isolated_mode = True if ('--no-build-isolation' in sys.argv) else None)

    if user_install:
        paths = (site.getusersitepackages(),)
    else:
        py_version = '%s.%s' % (sys.version_info[0], sys.version_info[1])
        paths = (s % (py_version) for s in (
            sys.prefix + '/lib/python%s/dist-packages',
            sys.prefix + '/lib/python%s/site-packages',
            sys.prefix + '/local/lib/python%s/dist-packages',
            sys.prefix + '/local/lib/python%s/site-packages',
            '/Library/Python/%s/site-packages',
        ))

    for path in paths:
        if os.path.exists(path):
            return path
    print('no installation path found', file=sys.stderr)
    return None


# Get the root folder of the CLASS installation -- this setup.py should be in that folder
root_folder = os.path.dirname(os.path.abspath(__file__))

# Recover the gcc compiler
GCCPATH_STRING = sbp.Popen(
    ['gcc', '-print-libgcc-file-name'],
    stdout=sbp.PIPE).communicate()[0]
GCCPATH = os.path.normpath(os.path.dirname(GCCPATH_STRING)).decode()

liblist = ["class"]
MVEC_STRING = sbp.Popen(
    ['gcc', '-lmvec'],
    stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec","m"]

# define absolute paths
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
            VERSION = line.split()[-1][2:-1]+".0" # Always set a sub-version number here for subsequent uploads (!)
            break

# Define cython extension and fix Python version
classy_ext = Extension("classy", [os.path.join("python", "classy.pyx")],
                       include_dirs=[np.get_include(), include_folder, heat_folder, recfast_folder, hyrec_folder],
                       libraries=liblist,
                       library_dirs=[root_folder, GCCPATH],
                       language="c++",
                       extra_compile_args=["-std=c++11"],
                       depends=["libclass.a","python/cclassy.pxd"]
                       )

classy_ext.cython_directives = {'language_level': "3" if sys.version_info.major>=3 else "2"}


# We need to put the actual files somewhere (e.g. the BBN file -- the easiest is just to copy the full class folder (although in the future we could imagine just copying the necessary files)

# 1. Check where to install class to
path_install = os.path.join(binaries_directory(),"class")
print("Selected installation path : ", path_install)


# 2. Check what files to include
def package_files(directory):
    paths = []
    direcs = []
    wanted_paths = {os.path.join(directory, d) for d in ["tools", "source", "main", "python", "include"]}
    for (path, directories, filenames) in os.walk(directory):
        # Only include those directories that we actually want
        if (path in wanted_paths or
            (path.startswith(os.path.join(directory,"external")) and not 'RealSpaceInterface' in path)):
          print("INCLUDING", path, filenames)
          for filename in filenames:
              paths.append(os.path.join(path, filename))
    return paths

pck_files = package_files(".")
pck_files.append("./Makefile")
# Debug print (only occurs when executed with the -v option)
print("Included files : ", pck_files)


# Make a custom builder in order to compile the C code as well, using the makefile
class classy_builder(build_ext):

    def build_extension(self, ext):
      # Make sure to put the current python version into the 'PYTHON' variable
      env = os.environ.copy()
      if not 'PYTHON' in env:
        env['PYTHON'] = sys.executable

      env['CLASSDIR'] = path_install

      # Compile the C code only
      returncode = sbp.call(["make","libclass.a","-j"], env=env)
      if returncode!=0:
        raise RuntimeError("Unknown error occurred -- the Makefile compilation of class failed (return code %i). Run the installation with '-v' and check why this is the case from the makefile command output"%returncode)
      super().build_extension(ext)

# Finally, perform the actual setup
setup(
    name='classy',
    version=VERSION,
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    cmdclass={'build_ext': classy_builder},
    ext_modules=[classy_ext],
    packages = ["class"],
    package_dir={"class":""},
    package_data={'class': pck_files},
    include_package_data=True
)
