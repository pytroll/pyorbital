.. pyorbital documentation master file, created by
   sphinx-quickstart on Mon Oct  3 08:48:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyorbital
=========

Pyorbital is a python package to compute orbital parameters for satellites from
TLE files as well as astronomical parameters of interest for satellite remote sensing.
Currently pyorbital only supports low earth orbit satellites.

TLE files
---------
Pyorbital has a module for parsing NORAD TLE-files

    >>> from pyorbital import tlefile
    >>> tle = tlefile.read('noaa 18', '/path/to/my/tle_file.txt')
    >>> tle.inclination
    99.043499999999995

If no path is given pyorbital tries to read the earth observation TLE-files from celestrak.com
    
Computing satellite postion
---------------------------
The orbital module enables computation of satellite position and velocity at a specific time:

    >>> from pyorbital.orbital import Orbital
    >>> from datetime import datetime
    >>> orb = Orbital("noaa 18")
    >>> now = datetime.utcnow()
    >>> # Get normalized position and velocity of the satellite:
    >>> orb.get_position(now)
    ([0.57529384846822862, 0.77384005228105424, 0.59301408257897559],
    [0.031846489698768146, 0.021287993461926374, -0.05854106186659274])
    >>> # Get longitude, latitude and altitude of the satellite:
    >>> orb.get_lonlatalt(now)
    (-1.1625895579622014, 0.55402132517640568, 847.89381184656702)

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

