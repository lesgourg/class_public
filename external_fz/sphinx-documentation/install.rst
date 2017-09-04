Installation Guide
==================

Prerequisites
-------------

For the use of the DarkAges package, you need a clean installation of Python_ (version 2.7)
with the the numpy_ and the scipy_ module.

.. warning:: Even though the code being written to be compatible with version 3.0,
   the compatibility is not fully tested yet

In addition to this minimal setting you will need to have the following two
packages to be installed

- **PyYaml**: This is for reading and saving additional options to the calculations
- **dill**: One of the key-points of the package is, that structures which are often used
  are stored in files and read again in later session rather to be recalculated every time,
  to save time

To test for the presence of the modules **numpy**,  **scipy**,
**dill**, **yaml** (PyYaml) on your machine, you can type

.. code::

   >>> import numpy
   >>> import scipy
   >>> import dill
   >>> import yaml

If one of these steps fails, go to the corresponding websites, and
follow the instructions (if you have the privilege to have the root
password on your machine, an `apt-get install python-numpy`,
`python-scipy` and `cython` will do the trick. Otherwise, all these
packages can also be downloaded and installed locally, with the
command :code:`python setup.py install --user`).

Using the modules
-----------------

Since this module is a pure Python module there is in principle no further need to install the code.

You are now free to use the Classes and methods provided in the package for example

.. code::

   >>> import DarkAges #import the full package
   >>> from DarkAges.common import time_at_z as t_of_z # time in s for a given redshift

and include them in your custom :code:`.py`-script.
For further informations on the modules and the methods and classes provided see the :ref:`documentation`

Alternatively we provide a command line script (see :ref:`using_the_command_line_script`)

.. _Python: http://www.python.org/
.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
