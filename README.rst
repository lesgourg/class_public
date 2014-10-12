==============================================
CLASS: Cosmic Linear Anisotropy Solving System
==============================================

:Author: Julien Lesgourgues

with several major inputs from other people, especially Thomas Tram,
as well as Benjamin Audren, Simon Prunet, Jesus Torrado, Miguel
Zumalacarregui, etc.

For download and information, see http://class-code.net


Compiling CLASS and getting started
-----------------------------------

(the information below can also be found on the webpage, just below
the download button)

After downloading the code, unpack the archive (tar -zxvf
class_v*.tar.gz), go to the class directory (cd class_v*/) and compile
(make clean; make class). If the first compilation attempt fails, you
may need to open the Makefile and adapt the name of the compiler
(default: gcc), of the optization flag (default: -O4) and of the
OpenMP flag (default: -fopenmp; this flag is facultative, you are free
to compile without OpenMP if you don't want parallel execution; note
that you need the version 4.2 or higher of gcc to be able to compile
with -fopenmp. Several details on the CLASS compilation are given on
the wiki page

https://github.com/lesgourg/class_public/wiki/Installation

(in particular, for compiling on Mac 10.9 Mavericks).

To check that the code runs, type:

    ./class explanatory.ini

The explanatory.ini file is a reference input file, containing and
explaning the use of all possible input parameters. We recommend to
read it, to keep it unchanged (for future reference), and to create
for your own purposes some shorter input files, containing only the
input lines which are useful for you. Input files must have a *.ini
extension.

If you want to play with the precision/speed of the code, you can use
one of the provided precision files (e.g. cl_permille.pre) or modify
one of them, and run with two input files, for instance:

    ./class test.ini cl_permille.pre

A simplified documentation can be found in the paper `CLASS I:
Overview <http://arxiv.org/abs/1104.2932>`_. On top of that, if you
wish to modify the code, you will find lots of comments directly into
the files. Other CLASS papers dedicated to various aspects of the code
are listed in the CLASS web page. Slides from CLASS-dedicated courses
can be seen at

http://lesgourg.web.cern.ch/lesgourg/class-tour/class-tour.html

To use CLASS from python, or ipython notebooks, or from the Monte
Python parameter extraction code, you need to compile not only the
code, but also its python wrapper. This can be done by typing just
'make' instead of 'make class'. More details on the wrapper and its
compilation are found on the wiki page

https://github.com/lesgourg/class_public/wiki

Plotting utility
----------------

Since version 2.3, the package includes an improved plotting script
called CPU.py (Class Plotting Utility), written by Benjamin Audren and
Jesus Torrado. It can plot the Cl's, the P(k) or any other CLASS
puput, for one or several models, as well as their ratio or percentage
difference. The syntax and list of available options is obtained by
typing 'pyhton CPU.py --help'. There is a similar script for MATLAB,
written by Thomas Tram. To use it, once in MATLAB, type 'help
plot_CLASS_output.m'

Developping the code
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
at least the paper `CLASS II: Approximation schemes
<http://arxiv.org/abs/1104.2933>`_. Feel free to cite more CLASS
papers!

Support
-------

To get support, please open a new issue on the

https://github.com/lesgourg/class_public

webpage!
