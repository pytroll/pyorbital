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
    >>> import datetime as dt
    >>> # Use current TLEs from the internet:
    >>> orb = Orbital("Suomi NPP")
    >>> now = dt.datetime.now(dt.timezone.utc)
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
    >>> import datetime as dt
    >>> orb = Orbital("Suomi NPP")
    >>> dtobj = dt.datetime(2015,2,7,3,0)
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



Instrument geolocation conventions
----------------------------------

The geolocation code builds a local orbital frame from the satellite position
and velocity, and rotates the instrument line of sight within it. Two choices in
that construction affect the ground solution, and both are selectable, because
changing them moves the geolocation of *existing* products.

**Nadir convention** -- which direction counts as "down":

``"legacy"`` (current default)
   The direction released versions of Pyorbital have always used: the
   normalised position of the ellipsoid subpoint of the antipode. The frame is
   left exactly as released, which means the nadir is *not* re-orthogonalised
   against the velocity.

``"geocentric"``
   Straight at the centre of the Earth. Validation of FY-3 MERSI against
   reference geolocation selects this one: the legacy convention leaves a
   latitude-dependent error going as ``sin(2 * latitude)``, about 0.3 km at 45
   degrees for an 850 km orbit and vanishing at the equator.

``"geodetic"``
   Along the ellipsoid normal. This departs furthest from the other two, by
   about 2.8 km at 45 degrees, and is *not* what the FY-3 validation supports.

**Rotation order** -- roll and pitch do not commute, and the cross-track field
of view carries the whole scan angle while pitch is a small attitude bias, so
the two orders diverge with scan angle: nothing at nadir, growing towards the
swath edge in proportion to the pitch bias.

``"legacy"`` (current default)
   Roll first, then pitch, as released. The along-track component of the
   pointing vector is compressed by ``cos(scan_angle)``.

``"pitch_first"``
   Pitch first, then roll, so the along-track component is independent of the
   scan angle. For AVHRR with a fitted pitch bias of 0.16 degrees the two
   differ by about 1.8 km at the swath edge; for FY-3F, whose bias is five
   times smaller, by about 0.37 km -- too little for that data to say which is
   right.

Selecting a convention
^^^^^^^^^^^^^^^^^^^^^^

Both can be passed directly to :func:`~pyorbital.geoloc.geolocate` and
:func:`~pyorbital.geoloc.compute_pixels`::

    >>> from pyorbital.geoloc import geolocate
    >>> lons, lats, alts = geolocate(orbital, scan_geometry, times,   # doctest: +SKIP
    ...                              nadir_convention="geocentric",
    ...                              rotation_order="pitch_first")

Libraries built on Pyorbital -- pygac and the readers using it, for instance --
call the geolocation without forwarding these arguments, so they can also be set
process-wide through the configuration, which reaches those callers with no
change to them::

    >>> import pyorbital                                              # doctest: +SKIP
    >>> with pyorbital.config.set(nadir_convention="geocentric",      # doctest: +SKIP
    ...                           rotation_order="pitch_first"):
    ...     reprocess()

or for a whole processing run, through the environment::

    PYORBITAL_NADIR_CONVENTION=geocentric PYORBITAL_ROTATION_ORDER=pitch_first ...

An explicit argument takes precedence over the configuration.

Reproducing existing products
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Leaving both unset keeps the released behaviour, so results do not move
underneath existing processing chains. Because that behaviour is measurably
less accurate, relying on the default raises a :exc:`DeprecationWarning`;
requesting ``"legacy"`` explicitly, or setting it in the configuration, is
treated as a deliberate choice and is silent.

Pinning both to ``"legacy"`` reproduces released Pyorbital to within floating
point (4e-11 degrees, well below a micrometre on the ground), which is what an
archive reprocessing needs::

    >>> with pyorbital.config.set(nadir_convention="legacy",          # doctest: +SKIP
    ...                           rotation_order="legacy"):
    ...     reprocess_the_archive()

Note this is a property of the pair. The legacy nadir also implies the released,
non-orthogonalised frame; mixing the legacy nadir with ``"pitch_first"`` does
*not* reproduce released output, and differs from it by kilometres once a pitch
bias is present.

Attitude biases fitted under one set of conventions do not carry over to
another: the rotation order in particular interacts with the pitch bias, so a
processing chain that estimates roll, pitch and yaw from ground control points
has to re-fit them if it changes conventions.

The defaults will change to the corrected conventions in a future release.


Computing astronomical parameters
---------------------------------
The astronomy module enables computation of certain parameters of interest for satellite remote sensing for instance the Sun-zenith angle:

    >>> from pyorbital import astronomy
    >>> import datetime as dt
    >>> utc_time = dt.datetime(2012, 5, 15, 15, 45)
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
