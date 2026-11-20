CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================

Authors: Julien Lesgourgues, Thomas Tram, Nils Schoeneberg

with several major inputs from other people, especially Benjamin
Audren, Simon Prunet, Jesus Torrado, Miguel Zumalacarregui, Francesco
Montanari, Deanna Hooper, Samuel Brieden, Daniel Meinert, Matteo Lucca, etc.

For download and information, see http://class-code.net


Compiling CLASS and getting started
-----------------------------------

(the information below can also be found on the webpage, just below
the download button)

Download the code from the webpage and unpack the archive (tar -zxvf
class_vx.y.z.tar.gz), or clone it from
https://github.com/lesgourg/class_public. Go to the class directory
(cd class/ or class_public/ or class_vx.y.z/) and compile (make clean;
make class). You can usually speed up compilation with the option -j:
make -j class.

You need a C compiler, a C++11 compiler (part of CLASS is C-style code
compiled as C++), and make. On macOS, install the Xcode Command Line
Tools (xcode-select --install); on Linux, install gcc/g++ and make from
your distribution. See "Choosing the compilers" below for using other
compilers, including those of a conda environment.

If the first compilation attempt fails, you may need to open the Makefile
and adapt the optimization flag (default: -O3) and the threading flag
(default: -pthread; CLASS no longer uses OpenMP by default, the line
with -fopenmp in the Makefile is commented out). Many more details on
the CLASS compilation are given on the wiki page

https://github.com/lesgourg/class_public/wiki/Installation

Choosing the compilers: the Makefile takes the C compiler from the
variable CC and the C++ compiler from the variable CXX. If you do not
set them, it uses "cc" and "c++", i.e. the system compilers (Apple clang on macOS, usually gcc/g++
on Linux). To use other compilers, set both variables, on the make
command line or in the environment (the environment also works for
"pip install", where the make command line is not accessible):

    make -j class CC=gcc-14 CXX=g++-14       # e.g. Homebrew gcc on macOS
    CC=gcc-14 CXX=g++-14 pip install .       # same, for the python wrapper

Always set CC and CXX together and from the same compiler family (clang
and clang++, or gcc and g++): the C++ standard library used by the C and
C++ parts of CLASS and by the python wrapper must be the same. After
changing compilers, run "make clean" first, because make does not
notice that the compiler changed.

In a conda environment, install the compilers and python packages with

    conda install compilers python numpy scipy cython setuptools wheel pip

(or "clang_osx-arm64 clangxx_osx-arm64", "clang_osx-64 clangxx_osx-64"
instead of "compilers" on macOS). Depending on the conda package,
activating the environment either defines CC and CXX (the
clang_osx-* packages), which the Makefile then uses, or defines nothing
but provides "cc" and "c++" (the "compilers" package, version 2.0), which
the Makefile uses by default. Either way the conda compilers are used,
with no further action. (Both cases were tested on macOS arm64. The
Fortran compiler of the "compilers" package is not needed by CLASS.)
BUILD_SYSTEM_FIXES.md explains these choices in detail, including what
was and was not tested.

To check that the code runs, type:

    ./class explanatory.ini

The explanatory.ini file is THE reference input file, containing and
explaining the use of all possible input parameters. We recommend to
read it, to keep it unchanged (for future reference), and to create
for your own purposes some shorter input files, containing only the
input lines which are useful for you. Input files must have a *.ini
extension. We provide an example of an input file containing a
selection of the most used parameters, default.ini, that you may use as a
starting point.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

The files *.pre are suppposed to specify the precision parameters for
which you don't want to keep default values. If you find it more
convenient, you can pass these precision parameter values in your *.ini
file instead of an additional *.pre file.

The automatically-generated documentation is located in

    doc/manual/html/index.html
    doc/manual/CLASS_manual.pdf

On top of that, if you wish to modify the code, you will find lots of
comments directly in the files.

Python
------

To use CLASS from python, or ipython notebooks, or from the Monte
Python parameter extraction code, you need to compile not only the
code, but also its python wrapper. This can be done by typing just
'make' instead of 'make class' (or for speeding up: 'make -j'), or with
'pip install .' (which also works as 'pip install classy' from PyPI).
The wrapper needs numpy, scipy and cython; in a conda environment,
install them first and use 'pip install --no-build-isolation .' so that
the build uses the environment's own numpy and cython. The compilers are
chosen as described in "Choosing the compilers" above. More
details on the wrapper and its compilation are found on the wiki page

https://github.com/lesgourg/class_public/wiki

Plotting utility
----------------

Since version 2.3, the package includes an improved plotting script
called CPU.py (Class Plotting Utility), written by Benjamin Audren and
Jesus Torrado. It can plot the Cl's, the P(k) or any other CLASS
output, for one or several models, as well as their ratio or percentage
difference. The syntax and list of available options is obtained by
typing 'pyhton CPU.py -h'. There is a similar script for MATLAB,
written by Thomas Tram. To use it, once in MATLAB, type 'help
plot_CLASS_output.m'

Developing the code
--------------------

If you want to develop the code, we suggest that you download it from
the github webpage

https://github.com/lesgourg/class_public

rather than from class-code.net. Then you will enjoy all the feature
of git repositories. You can even develop your own branch and get it
merged to the public distribution. For related instructions, check

https://github.com/lesgourg/class_public/wiki/Public-Contributing

Using the code
--------------

You can use CLASS freely, provided that in your publications, you cite
at least the paper `CLASS II: Approximation schemes <http://arxiv.org/abs/1104.2933>`. Feel free to cite more CLASS papers!

Support
-------

To get support, please open a new issue on the

https://github.com/lesgourg/class_public

webpage!
