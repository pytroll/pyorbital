.. pyorbital documentation master file, created by
   sphinx-quickstart on Mon Oct  3 08:48:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyorbital
=========

Pyorbital is a python package to compute orbital parameters for satellites from
TLE files as well as astronomical parameters of interest for satellite remote sensing.
Currently pyorbital only supports low earth orbit satellites.

Installation
------------
Pyorbital comes with a file platforms.txt that maps satellite name to NORAD identifier.
This file needs to be copied to the appropriate Satpy `etc` directory ($PPP_CONFIG_DIR).
It is wise to check it contains your satellites of interest. The NORAD identifier can
be found as the first number of each line in the Two-Line Elements (eg. from celestrak).

TLE files
---------
Pyorbital has a module for parsing NORAD TLE-files

    >>> from pyorbital import tlefile
    >>> tle = tlefile.read('noaa 18', '/path/to/my/tle_file.txt')
    >>> tle.inclination
    99.043499999999995

If no path is given pyorbital tries to read the earth observation TLE-files from celestrak.com


Content
=======

.. toctree::
   :maxdepth: 2

   sat_position
   astronomy
   snos
   api
                 
Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

