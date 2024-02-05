Pyorbital
=========

Pyorbital is a python package to compute orbital parameters for satellites from
TLE files as well as astronomical parameters of interest for satellite remote sensing.
Currently Pyorbital only supports low earth orbit satellites.


Installation
------------

Pyorbital is available from the Python Package Index (PyPI) via pip or from
the conda-forge conda channel. To install from PyPI in an existing environment:

.. code-block:: bash

   pip install pyorbital
   
Or in an existing conda-based environment:

.. code-block:: bash

   conda install -c conda-forge pyorbital

From Source
^^^^^^^^^^^

Pyorbital can also be installed from source. If you want to install pyorbital
from the latest in-development version on GitHub you can run:

.. code-block:: bash

   pip install git+https://github.com/pytroll/pyorbital.git
    
However, if you instead want to edit the source code and see the changes reflected
when you run the code you can clone the git repository and install it in
"editable" mode:

.. code-block:: bash

   git clone git://github.com/pytroll/pyorbital.git
   cd pyorbital
   pip install -e .


Add platform missing information
--------------------------------

Pyorbital comes with a file *platforms.txt* that maps a satellite name to the NORAD identifier.

This file already contain many low earth orbiting environmental or
meteorological satellites and thus likely be sufficient for your purpose.

But should it not contain your satellites of interest make a copy of the
`platforms.txt <https://github.com/pytroll/pyorbital/blob/main/pyorbital/etc/platforms.txt>`_
file and add the missing satellites and their NORAD identifiers and place
the file in the directory pointed to by :envvar:`PYORBITAL_CONFIG_PATH`.

The NORAD identifier can be found as the first number of each line in the
Two-Line Elements files (eg. from `celestrak`_).

Pyorbital comes with a small script ``check_platform.py`` to check whether a
satellite is already supported.

.. code::

   python -m pyorbital.check_platform -s NOAA-21

   [INFO: 2023-01-22 21:20:25 : pyorbital.tlefile] Satellite NOAA-21 is supported. NORAD number: 54234
   [INFO: 2023-01-22 21:20:25 : pyorbital.tlefile] Satellite names and NORAD numbers are defined in /path/to/pyorbital/etc/directory/platforms.txt


TLE files
---------
Pyorbital has a module for parsing NORAD TLE-files

    >>> from pyorbital import tlefile
    >>> tle = tlefile.read('noaa 18', '/path/to/my/tle_file.txt')
    >>> tle.inclination
    99.043499999999995

If no path is provided pyorbital first tries to read any local TLE files defined by the
environment variable :envvar:`TLES` giving a glob pattern that can be used to retrieve all relevant files:

.. code::

   TLES=/path/to/tle_files/*/tle*txt

If this variable is not set Pyorbital will try get the earth observation TLE files over the internet
from `celestrak`_. Note this downloading only happens if no
specific TLE file is provided or if the :envvar:`TLES` environment variable is not set.


TLE download and database
^^^^^^^^^^^^^^^^^^^^^^^^^

The historical TLE files can be requested from
`celestrak's request page <https://celestrak.com/NORAD/archives/request.php>`_.

There is also a script, ``fetch_tles.py``, that can be used to collect
TLE data from several locations. The currently supported locations
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

    >>> snpp = Orbital('Suomi NPP', tle_file='/path/to/tle/files/tle-20150207.txt')
    >>> snpp.get_lonlatalt(dtobj)
    (105.37373804512762, 79.160752404540133, 838.94605490133154)

If we take a TLE from one week earlier we get a slightly different result:

    >>> snpp = Orbital('Suomi NPP', tle_file='/path/to/tle/files/tle-20150131.txt')
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


.. envvar:: PYORBITAL_CONFIG_PATH

   It is possible (but not mandatory) to define this environment variable to
   have full control of certain static data used by Pyorbital:

   Pyorbital comes with a file *platforms.txt* that maps a satellite name to the
   NORAD identifier. This internal file is accessed by Pyorbital without the
   user having to do anything. But if you need to change or update this file
   you can make your own copy and place in the directory pointed to by this
   environment variable.

.. envvar:: TLES

   Two Line Element (TLE) files are accessed automatically over the internet
   without the user having to do anything. When doing that Pyorbital will fetch
   the most recent TLE data which may not be the most optimal for historic data
   for instance. Also, it may not be sustainable in a production environment.

   However, it is possible to let Pyorbital look for the necessary and more
   optimal TLE data locally, by specifying locations where such local TLE
   files are located. If the TLES environment variable is set to a glob pattern to
   local locations, Pyorbital will first search for the needed TLEs
   there. This can both be useful in an operational setup where access to the
   internet is restricted, and when processing old/historic satellite data.

   It is possible (but not mandatory) to define this environment variable.


API
---

Orbital computations
^^^^^^^^^^^^^^^^^^^^

.. automodule:: pyorbital.orbital
   :members:
   :undoc-members:

TLE handling
^^^^^^^^^^^^

.. automodule:: pyorbital.tlefile
   :members:
   :undoc-members:

Astronomical computations
^^^^^^^^^^^^^^^^^^^^^^^^^

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



.. _celestrak: Celestrak <https://celestrak.com>
.. _github: http://github.com/pytroll/pyorbital
