from setuptools import setup, Extension
from Cython.Distutils import build_ext

import numpy as np
import os
import subprocess as sbp
import sys

def binaries_directory():
    path = "."

    if '--target' in sys.argv:
      path = os.path.abspath(sys.argv[sys.argv.index('--target')+1])

    if 'editable_wheel' in sys.argv:
      path = os.path.dirname(os.path.abspath(__file__))

    if path == ".":
      print('no installation path found', file=sys.stderr)
    return path

bin_dir = binaries_directory()
print("BIN_DIR",bin_dir)

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
classy_ext = Extension("classy._classy", [os.path.join("python", "classy.pyx")],
                       include_dirs=[np.get_include(), include_folder, heat_folder, recfast_folder, hyrec_folder],
                       libraries=liblist,
                       library_dirs=[root_folder, GCCPATH],
                       language="c++",
                       extra_compile_args=["-std=c++11"],
                       depends=["libclass.a","python/cclassy.pxd"]
                       )

classy_ext.cython_directives = {'language_level': "3" if sys.version_info.major>=3 else "2"}




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

      # We need to put the actual files somewhere (e.g. the BBN file -- the easiest is just to copy the full class folder (although in the future we could imagine just copying the necessary files)

      # 1. Check where to install class to
      print("Running installation with arguments = ",sys.argv)
      path_install = binaries_directory()
      print("Selected corresponding installation path : ", path_install)

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
    packages = ["classy"],
    package_dir={"classy":"."},
    package_data={'classy': pck_files},
    include_package_data=True,
    zip_safe=False
)
