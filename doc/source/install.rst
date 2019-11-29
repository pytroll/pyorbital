=========================
Installation Instructions
=========================

Note
----
Pyorbital comes with a file platforms.txt that maps satellite name to NORAD identifier.
This file needs to be copied to the appropriate Satpy `etc` directory ($PPP_CONFIG_DIR).
It is wise to check it contains your satellites of interest. The NORAD identifier can
be found as the first number of each line in the Two-Line Elements (eg. from celestrak).

Pip-based Installation
======================

Pyorbital is available from the Python Packaging Index (PyPI). A sandbox
environment for `pyorbital` can be created using
`Virtualenv <http://pypi.python.org/pypi/virtualenv>`_.

To install the `pyorbital` package and the python dependencies:

.. code-block:: bash

    $ pip install pyorbital


Conda-based Installation
========================

Starting with version 1.3.1, Pyorbital is available from the conda-forge channel. If
you have not configured your conda environment to search conda-forge already
then do:

.. code-block:: bash

    $ conda config --add channels conda-forge

Then to install Pyorbital in to your current environment run:

.. code-block:: bash

    $ conda install pyorbital
