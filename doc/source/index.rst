.. pyorbital documentation master file, created by
   sphinx-quickstart on Mon Oct  3 08:48:29 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyorbital
=========

Pyorbital is a python package to compute orbital parameters for satellites from
TLE files.

Usage
-----

Here's a typical usage for the package.

  >>> from pyorbital.orbital import Orbital
  >>> from datetime import datetime
  >>> orb = Orbital("noaa 18")
  >>> now = datetime.utcnow()
  >>> # Get normalized position of the satellite:
  >>> orb.get_position(now)
  ([0.57529384846822862, 0.77384005228105424, 0.59301408257897559],
  [0.031846489698768146, 0.021287993461926374, -0.05854106186659274])
  >>> # Get longitude, latitude and altitude of the satellite:
  >>> orb.get_lonlatalt(now)
  (-1.1625895579622014, 0.55402132517640568, 847.89381184656702)

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

