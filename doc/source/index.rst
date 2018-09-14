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
This file needs to be copied to the appropriate satpy etc directory ($PPP_CONFIG_DIR).
It is wise to check it contains your satellites of interest. The NORAD identifier can
be found as the first number of each line in the Two-Line Elements (eg. from celestrak).

TLE files
---------
Pyorbital has a module for parsing NORAD TLE-files


.. _github: http://github.com/pytroll/pyorbital

    >>> from pyorbital.orbital import Orbital
    >>> from datetime import datetime
    >>> # Use current TLEs from the internet:
    >>> orb = Orbital("Suomi NPP")
    >>> now = datetime.utcnow()
    >>> # Get normalized position and velocity of the satellite:
    >>> orb.get_position(now)
    (array([-0.20015267,  0.09001458,  1.10686756]),
     array([ 0.06148495,  0.03234914,  0.00846805]))
    >>> # Get longitude, latitude and altitude of the satellite:
    >>> orb.get_lonlatalt(now)
    (40.374855865574951, 78.849923885700363, 839.62504115338368)
 

Use actual TLEs to increase accuracy
------------------------------------

    >>> from pyorbital.orbital import Orbital
    >>> from datetime import datetime
    >>> orb = Orbital("Suomi NPP")
    >>> dtobj = datetime(2015,2,7,3,0)
    >>> orb.get_lonlatalt(dtobj)
    (152.11564698762811, 20.475251739329622, 829.37355785502211)

But since we are interesting knowing the position of the Suomi-NPP more than
two and half years from now (September 26, 2017) we can not rely on the current
TLEs, but rather need a TLE closer to the time of interest:

    >>> snpp = Orbital('Suomi NPP', tle_file='/data/lang/satellit/polar/orbital_elements/TLE/201502/tle-20150207.txt')
    >>> snpp.get_lonlatalt(dtobj)
    (105.37373804512762, 79.160752404540133, 838.94605490133154)

If we take a TLE from one week earlier we get a slightly different result:

    >>> snpp = Orbital('Suomi NPP', tle_file='/data/lang/satellit/polar/orbital_elements/TLE/201501/tle-20150131.txt')
    >>> snpp.get_lonlatalt(dtobj)
    (104.1539184988462, 79.328272480878141, 838.81555967963391)


>>>>>>> master

.. toctree::
   :maxdepth: 2

   tle
   satellite_position
   astronomy
   moon_calculations
   api
              

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

