Moon phase and position with Pyorbital
======================================

It is possibe to calculate the phase and position of the moon and any given
time in the past or future. This is done by an implementation of the
astronimical methods described at Stjaernhimlen_. The calculations have been
compared to the pyephem library, and deviations are small, see below.


Computing the phase of the moon
-------------------------------

  >>> from datetime import datetime
  >>> from pyorbital.moon_phase import moon_phase
  >>> time_t = datetime(2011, 12, 1, 12)
  >>> print moon_phase(time_t)
  0.409752921579

Computing the position of the moon
----------------------------------

  >>> from datetime import datetime
  >>> from pyorbital import planets
  >>> time_t = datetime(2016, 7, 30, 0, 0)
  >>> moon = planets.Moon(time_t)
  >>> lon = 20.0
  >>> lat = 65.0
  >>> rasc, decl, alt, azi = moon.topocentric_position(lon, lat)
  >>> print alt, azi
  6.33389454706 61.6795817556


Accuracy
--------

In lack of an absolute truth or reference the moon calclations with Pyorbital
have been compared to pyephem. There are indeed deviations, but for the moon
phase they are in general rather small. See image below, where we compare the
moon phase over 416 days starting from December 1st 2015. As seen from the figure the deviations are within 0.3 %.

  .. image:: _static/moonphase_compare.png

For moon height and azimuth differences are larger, as seen from the figures
below. We have calculated the position of the moon relative to the City of
Norrköping, Sweden, over four months from March 7, 2012. The deviations are
within 4 degrees for the height, and within 9 degrees for the azimuth.


  .. image:: _static/moonheight_compare.png
  .. image:: _static/moonazimuth_compare.png



.. _`Stjaernhimlen`:   http://www.stjarnhimlen.se/comp/ppcomp.html
