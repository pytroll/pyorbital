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

    >>> from pyorbital import tlefile
    >>> tle = tlefile.read('noaa 18', '/path/to/my/tle_file.txt')
    >>> tle.inclination
    99.043499999999995

If no path is given pyorbital tries to read the earth observation TLE-files from celestrak.com

TLE download and database
~~~~~~~~~~~~~~~~~~~~~~~~~

The historical TLE files can be requested from
`celestrak <https://celestrak.com/NORAD/archives/request.php>`_.

There is also a script, ``fetch_tles.py``, that can be used to collect
TLE data from several locations.  Then currently supported locaions
are:

* generic network locations without login
* Space-Track (login credentials needed)
* local files

The data are saved in a SQLite3 database, and can be written to a file
after each run.  To see configuration options, see the example
configuration in ``examples/tle.yaml``.

Computing satellite position
----------------------------
The orbital module enables computation of satellite position and velocity at a specific time:

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

But since we are interested in knowing the position of the Suomi-NPP more than
two and half years from now (September 26, 2017) we can not rely on the current
TLEs, but rather need a TLE closer to the time of interest:

    >>> snpp = Orbital('Suomi NPP', tle_file='/data/lang/satellit/polar/orbital_elements/TLE/201502/tle-20150207.txt')
    >>> snpp.get_lonlatalt(dtobj)
    (105.37373804512762, 79.160752404540133, 838.94605490133154)

If we take a TLE from one week earlier we get a slightly different result:

    >>> snpp = Orbital('Suomi NPP', tle_file='/data/lang/satellit/polar/orbital_elements/TLE/201501/tle-20150131.txt')
    >>> snpp.get_lonlatalt(dtobj)
    (104.1539184988462, 79.328272480878141, 838.81555967963391)



Computing astronomical parameters
---------------------------------
The astronomy module enables computation of certain parameters of interest for satellite remote sensing for instance the Sun-zenith angle:

    >>> from pyorbital import astronomy
    >>> from datetime import datetime
    >>> utc_time = datetime(2012, 5, 15, 15, 45)
    >>> lon, lat = 12, 56
    >>> astronomy.sun_zenith_angle(utc_time, lon, lat)
    62.685986438071602

API
---

Orbital computations
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pyorbital.orbital
   :members:
   :undoc-members:

TLE handling
~~~~~~~~~~~~

.. automodule:: pyorbital.tlefile
   :members:
   :undoc-members:

Astronomical computations
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: pyorbital.astronomy
   :members:
   :undoc-members:


.. Contents:
   .. toctree::
      :maxdepth: 2
   Indices and tables
   ==================
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

