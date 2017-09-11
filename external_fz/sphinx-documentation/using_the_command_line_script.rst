.. _using_the_command_line_script:

Working with the command-line script
====================================

Since every working session with the :code:`DarkAges`-package has
basically the same structure, namely

#. Setting up the spectra of electrons/positrons, photons, and inefficcient
   particles for the given model parameters, by reading them from a table and
   interpolating and/or calculate them analytically.
#. Assign a redshift dependence on the spectra according to the 
   energy injection history in question and define an instance of the
   :class:`model <DarkAges.model.model>`-class.
#. Calculating the *effective efficiency factor* :math:`f(z)` for the given
   :class:`model <DarkAges.model.model>` with a given set of
   transfer-functions :math:`T_{klm}` provided in an initialized instance
   of the :class:`transfer <DarkAges.transfer.transfer>`-class.
#. Store the resulting table in a file or print it on a screen.

we provide a command-line scripts which does most the steps automatically. It is 
located under :code:`./bin/DarkAges` in the root-directory of the package. Even though
its location may propose that it is a binary file, it is indeed an executable
python script, saying that it can be executed in both ways, like

.. code::

	$ ./bin/DarkAges ...

or directly in the manner of an python script with

.. code::

	$ python ./bin/DarkAges ...

Especially if you are using more than one installations of python we advise you to 
execute the script in the latter way, with repalcing :code:`python` by the path to
the installation of python you want to use, if needed.

.. note::

   Please ensure that the script has the execution flag enabled, if you want to execute the script
   directly. If it is not set, run

   .. code::

      $ chmod +x ./bin/DarkAges

The structure command-line script is in principle that it parses the input parameters
given in the command-line, runs some basic consistency checks, and depending on the 
input-values executes a given routine (this routines are part of the 
:mod:`recipes <DarkAges.recipes>`-module and can also be used in a custom
python script using the DarkAges-package).  

Basic parameters
----------------

This parameters specify the details of the DM-model in question, the injection history, and the
cosmological background parameters

:code:`--mass 65`
   The mass of the DM candidate (*Here: 65*) in units of :math:`\mathrm{GeV}`, :math:`\mathrm{g}`, :math:`\mathrm{M_\mathrm{sun}}` depending 
   on the specified type injection history. This parameter is obligatory for the *"Live calculation mode"* (see below).  

:code:`--tdec 1e17`
   Lifetime of the DM candidate in units of :math:`\mathrm{s}^{-1}`.
   This is obligatory for the :code:`decay`-history. For the other cases this parameter will not be considered.

:code:`--hist annihilation`
   The injection history. This is needed to apply the correct redshift dependence on the injection spectra and
   the scaling of the number density with redshift. This parameter is optional. Per default :code:`annihilation`
   is taken. The valid options are :code:`annihilation`, :code:`decay`, or :code:`PBH`

:code:`--use-background 67 0.3 8e-5`
   Throughout the code, for example for the convolution with the transfer functions,
   the value of :code:`H(z)` is needed, which depends on the values of :math:`H_0`,
   :math:`\Omega_\mathrm{mat.}`, and :math:`\Omega_\mathrm{rad.}`. To perform a consistent analyses with CLASS this
   values can be passed to the calculations. The value of :math:`H_0` needs to be given in units of :math:`\frac{\mathrm{km}}{\mathrm{Mpc}\,\mathrm{s}}`.

:code:`--extra-options options.yaml`
   File with additional parameters and options, like precision parameters, passed to the methods of :code:`DarkAges` (as part of the
   :code:`DarkOptions`-structure). This file needs to written in the structure of YAML.


The three execution modes
-------------------------
 
The command-line script has basically three different modes in which it can be executed.

The *"Live calculation mode"*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this mode all the steps given above are done in one session. It is basically started with
