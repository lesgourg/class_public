crash-less & (partial) _h_-less CLASS {#mainpage}
=============================================
(Scroll down for original CLASS v3.2.0 README.md contents)

Author: Kevin Croker

This verison forked from v3.2.0.
The CLASS team did a bunch of reworks, like improving the Python side of the build chain and
removing OpenMP in subsequent versions, I've not ported the features below to that stuff yet.

Crash-less Explanation
------------------------
When CLASS (the C stuff) errors out because of a bad cosmology, it didn't free its allocated memory.
This was no problem if just called from the command line because the OS cleans up for you.
But if called from the Python wrapper, this became a severe memory leak.
If driven from MontePython for MCMC, available memory could be exhausted in seconds if the MCMC got into some squirelly place in cosmological parameter space.
This was not a rare bug, it was around for many years and people complained about it on github, and annoyingly worked around it by just restarting MontePython, as many times as necessary.

This version of CLASS uses a very lightweight doubly linked list memory tracker so that, when failure happens, all allocated memory can be freed reliably.
In my own use, MontePython can now run indefinitely and other diagnostics show stable memory usage over many hours of operation.
(if you want to lol, run valgrind on CLASS ;) )

_h_-less Explanation
------------------------
Under the hood, CLASS integrates physical densities in 1/Mpc<sup>2</sup> units.
Yet, for accounting, closing, and in various other places convenience, it works with big &Omega;'s and _H_<sub>0</sub>.
Although this is the standard way of "thinking" about Friedmann models, this is actually really bad for searching model space because it binds a present-day scale (_H_<sub>0</sub>) to the code that establishes initial conditions.
For example, by writing &Omega;<sub>b</sub> _H_<sub>0</sub><sup>2</sup>/_a_<sub>i</sub>^<sup>3</sup> to get the baryon density at _a_<sub>i</sub>, we've just projected backwards assuming nothing squirrely is happening across the intervening 13.8Gyr.
Its a bold strategy Cotton, lets see if it pays off for 'em.

Anyway.
In this version of CLASS, the `background` and `thermodynamics` modules have been modified to work without specification of an _H_<sub>0</sub> or an _h_.
Its an option, that can be flagged with
```
without_h = yes
```
What then happens is that you are only permitted to specify the little omegas (&omega;'s).
These are already physical quantities, because they have a density scale baked into them.
CLASS then interprets the &omega; quantities as densities that _can_ be projected backwards in time for initial conditions.
Anywhere in the code that made projection assumptions at intermediate scale factors, presumably for simplicity, has been adjusted to use the actual integrated densities determined by CLASS in the background module.

The motivation for this is that CMB experiments are telling you about early universe physics, before squirrely late-time things could be happening.
So the best-fit parameters should be describing the early-universe.
Other parameters can describe departures from the post-recombination conditions established by the little &omega;'s.
The upshot of all this is that early-time and late-time physical processes determine the resultant _H_<sub>0</sub>, so it obviously becomes a derived parameter, as it should be.

TLDR - if `without_h` is enabled, then the &omega;'s describe the universe that would result today if late-time physical processes do not alter the background expansion history.  Hubble can no longer be specified, as it becomes determined by the expansion history from initial conditions set by the &omega;'s projected backwards.

**CAVEAT** - Only implemented in `background` and `thermodynamics` for now.

(Below is the CLASS boilerplate)

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
make -j class. If the first compilation attempt fails, you may need to
open the Makefile and adapt the name of the compiler (default: gcc),
of the optimization flag (default: -O4 -ffast-math) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution; note that
you need the version 4.2 or higher of gcc to be able to compile with
-fopenmp). Many more details on the CLASS compilation are given on the
wiki page

https://github.com/lesgourg/class_public/wiki/Installation

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

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
'make' instead of 'make class' (or for speeding up: 'make -j'). More
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
