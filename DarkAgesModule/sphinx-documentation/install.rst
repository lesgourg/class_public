Installation Guide
==================

Prerequisites
-------------

For the use of the DarkAges package, you need a clean installation of Python_
(version 2.7 or version 3.x) with the the numpy_ and the scipy_ module.

Since version 1.1.0 the code is compatible to work with both common major
versions 2.7 and 3.x.
To ensure upward compatibility the **future** package is required.

In addition to this minimal setting you will need to have the following two
packages to be installed

- **PyYaml**: This is for reading and saving additional options to the calculations
- **dill**: One of the key-points of the package is, that structures which are often used
  are stored in files and read again in later session rather to be recalculated every time,
  to save time

All requirements (if not already satisfied) can be installed via
:code:`pip install -r requirements.txt`.
In case you do not have administrator rights on the machine, add :code:`--user`
to the command above to install the required packages locally.
Please make sure that you are using the correct version of :code:`pip` if you
have multiple installations of Python_ on your machine.

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
