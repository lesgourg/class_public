CLASS: Cosmic Linear Anisotropy Solving System  
==============================================

Author: Julien Lesgourgues

_This manual is under construction; this is only a provisional version. The definitive version will be made available soon, as well as all the necessary documentation to generate new versions of the manual. Currently the introduction is outdated and the definitions for some specific variables in the header files are missing. There are also some unresolved formatting issues in the documentation for harmonic.c and transfer.c, which will be corrected shortly._

Overall architecture of `class`
==========================================

The seven-module backbone
-------------------------

The purpose of `class` consists in computing some power
spectra for a given set of cosmological parameters. This task can be
decomposed in few steps or modules:

1.  compute the evolution of cosmological background quantitites.

2.  compute the evolution of thermodynamical quantitites (ionization
    fractions, etc.)

3.  compute the evolution of source functions \f$S(k,\eta)\f$ (by
    integrating over all perturbations).

4.  compute Bessel functions (in order to go from Fourier to harmonic
    space).

5.  compute transfer functions \f$\Delta_l(k)\f$ (unless one needs only
    Fourier spectra \f$P(k)\f$’s and no harmonic spectra \f$C_l\f$’s).

6.  compute the primordial spectrum for scalars, tensors, etc.
    (straightforward if the input consists in spectral parameters \f$A_s\f$,
    \f$n_s\f$, \f$r\f$, ..., but this module will incorporate the option of
    integrating over inflationary perturbations).

7.  compute power spectra \f$C_l\f$’s and/or \f$P(k)\f$’s.

In `class`, each of these steps is associated with a
structure:

1.  `struct background ` for cosmological background,

2.  `struct thermodynamics ` for thermodynamics,

3.  `struct perturbations ` for source functions,

4.  `struct bessels ` for bessel functions,

5.  `struct transfer ` for transfer functions,

6.  `struct primordial ` for primordial spectra,

7.  `struct harmonic ` for output spectra.

A given structure contains “everything concerning one step that the
subsequent steps need to know” (for instance, everything about source
functions that the transfer module needs to know). In particular, each
structure contains one array of tabulated values (background quantitites
as a function of time, thermodynamical quantitites as a function of
redshift, sources as a function of \f$(k, \eta)\f$, etc.). It also contains
information about the size of this array and the value of the index of
each physical quantity, so that the table can be easily read and
interpolated. Finally, it contains any derived quantity that other
modules might need to know. Hence, the comunication from one module A to
another module B consists in passing a pointer to the structure filled
by A, and nothing else.

Each structure is defined and filled in one of the following modules
(and precisely in the order below):

1.  `background.c `

2.  `thermodynamics.c `

3.  `perturbations.c `

4.  `bessel.c `

5.  `transfer.c `

6.  `primordial.c `

7.  `harmonic.c `

Each of these modules contains at least three functions:

-   *module*\_`init(...)`

-   *module*\_`free()`

-   *module*\_*something*\_`at`\_*somevalue*`(...)`

The first function allocates and fills each structure. This can be done
provided that the previous structures in the hierarchy have been already
allocated and filled. In summary, calling one of
*module*\_`init(...)` amounts in solving
entirely one of the steps 1 to 7.

The second function deallocates the fields of each structure. This can
be done optionally at the end of the code (or, when the code is embedded
in a sampler, this *must* be done between each execution of
`class`, and especially before calling
*module*\_`init(...)` again with different input
parameters).

The third function is able to interpolate the pre-computed tables. For
instance, `background\_init()` fills a table of background
quantitites for discrete values of conformal time \f$\eta\f$, but
`background\_at\_eta(eta, \* values)` will return these
values for any arbitrary \f$\eta\f$.

Note that functions of the type
*module*\_*something*\_`at`\_*somevalue*`(...)`
are the only ones which are called from another module, while functions
of the type *module*\_`init(...)` and
*module*\_`free()` are the only one called by
the main executable. All other functions are for internal use in each
module.


Input
-----

There are two types of input:

1.  “precision parameters” (controlling the precision of the output and
    the execution time),

2.  “input parameters” (cosmological parameters, flags telling to the
    code what it should compute, ...)

All “precision parameters” have been grouped in a single structure
`struct precision`. The code contains *no other
arbitrary numerical coefficient*. This structure is initialized
in a simple module `precision.c` by the function
`precision\_init()`. Nothing is allocated dynamically in this
function, so there is no need for a `precision\_free()`
function.

Each “input parameter” refers to one particular step in the computation:
background, thermodynamics, perturbations, etc. Hence they are defined
as part of the corresponding structure. Their values are assigned in a
simple module `input.c`, by a function
`input\_init(...)` which has a pointer towards each structure
in its list of arguments. Hence, when a given function
*module*\_`init(...)` is called, the
corresponding structure already contains input parameters; the function
fills the rest of this structure. The function
`input\_init(...)` does not allocate any field dynamically,
so there is no need for an `input\_free()` function.

Output
------

A simple module `output.c` writes the final results in files.
The name of the files are considered as input parameters making part of
a small structure `struct output`. Like for all other input
parameters, these names are assigned inside the function
`input\_init(...)`. Again this structure contains no
dynamically allocated quantitites, so there is no need for an
`output\_free()` function.

Summary
-------

We hope that after this short overview, it is clear for the reader that
the main executable of `class` should consist only in the
following lines (not including comments and error-management lines):

For a given purpose, somebody could only be interested in the
intermediate steps (only background quantitites, only the
thermodynamics, only the perturbations and sources, etc.) It is then
straightforward to truncate the full hierarchy of modules 1, ... 7 at
some arbitrary order. We provide several “reduced executables”
`test`\_*module* achieving precisely this.

Note also that if `class` is embedded in a parameter sampler
and only “fast” parameters are varied (i.e., parameters related to the
primordial spectra), then it is only necessary to repeat the following
steps after `output\_init(...)`:

`spectra\_free()

primordial\_free()

input\_init(&ba,&th,&pt,&bs,&tr,&pm,&hr,&op)

primordial\_init(&pt,&pr,&pm)

spectra\_init(&pt,&tr,&pm,&hr)

output\_init(&pt,&tr,&hr,&op)

`

General principles
==================

Flexibility
-----------

Explain allocation of indices, ...

Control of precision
--------------------

Explain precision structure, ...

Control of errors
-----------------
