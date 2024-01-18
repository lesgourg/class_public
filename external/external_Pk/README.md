The `external_Pk` mode
======================

* Author: Jesus Torrado (torradocacho [@] lorentz.leidenuniv.nl)
* Date:   2013-12-20


Introduction
------------

This mode allows for an arbitrary primordial spectrum `P(k)` to be calculated by an external command and passed to CLASS. That external command may be anything that can be run in the shell: a python script, some compiled C or Fortran code... This command is executed from within CLASS, and CLASS is able to pass it a number of parameters defining the spectrum (an amplitude, a tilt...). Those parameters can be used in a Markov chain search performed by MontePython.

This mode includes the simple case of a precomputed primordial spectrum stored in a text file. In that case, the `cat` shell command will do the trick (see below).

Currently, scalar and tensor spectra of perturbations of adiabatic modes are supported.


Use case #1: reading the spectrum from a table
----------------------------------------------

In this case, say the file with the table is called `spectrum.txt`, located under `/path/to`, simply include in the `.ini` file

    command = cat path/to/spectrum.txt
		
It is necessary that 1st 4 characters are exactly `cat `.


Use case #2: getting the spectrum from an external command
----------------------------------------------------------

Here an external command is called to generate the spectrum; it may be some compiled C or Fortran code, a python script... This command may be passed up to 10 floating point arguments, named `custom1` to `custom10`, which are assigned values inside the `.ini` file of CLASS. The `command` parameter would look like

    command = /path/to/example.py

if it starts with `#/usr/bin/python`, otherwise

    command = python /path/to/example.py

As an example of the 1st use case, one may use the included script `generate_Pk_example.py`, which implements a single-field slow-roll spectrum without running, and takes 3 arguments:
* `custom1` -- the pivot scale (`k_0 = 0.05 1/Mpc` for Planck).
* `custom2` -- the amplitude of the scalar power spectrum.
* `custom3` -- the scalar spectral index.

In order to use it, the following lines must be present in the parameter file:

    P_k_ini type = external_Pk
    command = /path/to/CLASS/external_Pk/generate_Pk_example.py
    custom1 = 0.05
    custom2 = 2.2e-9
    custom3 = 1.

Defined or not (in that case, 0-valued), parameters from `custom4` to `custom10` will be passed to the example script, which should ignore them. In this case, CLASS will run in the shell the command

    /path/to/CLASS/external_Pk/generate_Pk_example.py 0.05 2.2e-9 1. 0 0 0 0 0 0 0

If CLASS fails to run the command, try to do it directly yourself by hand, using exactly the same string that was given in `command`.


Output of the command / format of the table
-------------------------------------------

The command must generate an output separated into lines, each containing a tuple (`k`, `P(k)`). The following requirements must be fulfilled:

* Each line must contain 2 (3, if tensors) floating point numbers: `k` (in `1/Mpc` units) and `P_s(k)` (and `P_t(k)`, if tensors), separated by any number of spaces or tabs. The numbers can be in scientific notation, e.g. `1.4e-3`.

* The lines must be sorted in increasing values of `k`.

* There must be at least two points `(k, P(k))` before and after the interval of `k` requested by CLASS, in order not to introduce unnecessary interpolation error. Otherwise, an error will be raised. In most of the cases, generating the spectrum between `1e-6` and `1 1/Mpc` should be more than enough.


Precision
---------

This implementation properly handles double-precision floating point numbers (i.e. about 17 significant figures), both for the input parameters of the command and for the output of the command (or the table).

The sampling of `k` given by the command (or table) is preserved to be used internally by CLASS. It must be fine enough a sampling to clearly show the features of the spectrum. The best way to test this is to plot the output/table and check it with the naked eye.

Another thing to have in mind arises at the time of convolving with the transfer functions. Two precision parameters are implied: the sampling of `k` in the integral, given by `k_step_trans`, and the sampling of the transfer functions in `l`, given by `l_logstep` and `l_linstep`. In general, it will be enough to reduce the values of the first and the third parameters. A good start is to give them rather small values, say `k_step_trans=0.01` and `l_linstep=1`, and to increase them slowly until the point at which the effect of increasing them gets noticeable.


Parameter fit with MontePython
------------------------------

(MontePython)[http://montepython.net/] is able to interact with the `external_Pk` mode transparently, using the `custom` parameters in an MCMC fit. One must just add the appropriate lines to the input file of MontePython. For our example, if we wanted to fit the amplitude and spectral index of the primordial spectrum, it would be:

    data.cosmo_arguments['P_k_ini type'] = 'external_Pk'
    data.cosmo_arguments['command'] = '/path/to/CLASS/external_Pk/generate_Pk_example.py'
    data.cosmo_arguments['custom1'] = 0.05                                   # k_pivot
    data.parameters['custom2']      = [ 2.2,  0, -1,  0.055, 1.e-9, 'cosmo'] # A_s
    data.parameters['custom3']      = [  1.,  0, -1, 0.0074,     1, 'cosmo'] # n_s

Notice that since in our case `custom1` represents the pivot scale, it is passed as a (non-varying) argument, instead of as a (varying) parameter.

In this case, one would not include the corresponding lines for the primordial parameters of CLASS: `k_pivot`, `A_s`, `n_s`, `alpha_s`, etc. They would simply be ignored.


Limitations
-----------

* So far, this mode cannot handle vector perturbations, nor isocurvature initial conditions.
* The external script knows nothing about the rest of the CLASS parameters, so if it needs, e.g., `k_pivot`, it should be either hard coded, or its value passed as one of the `custom` parameters.
