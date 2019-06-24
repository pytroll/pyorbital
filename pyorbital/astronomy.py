#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2013

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Astronomy module.
Parts taken from http://www.geoastro.de/elevaz/basics/index.htm
"""

import numpy as np

from pyorbital import dt2np

F = 1 / 298.257223563  # Earth flattening WGS-84
A = 6378.137  # WGS84 Equatorial radius
MFACTOR = 7.292115E-5


def jdays2000(utc_time):
    """Get the days since year 2000.
    """
    return _days(dt2np(utc_time) - np.datetime64('2000-01-01T12:00'))


def jdays(utc_time):
    """Get the julian day of *utc_time*.
    """
    return jdays2000(utc_time) + 2451545


def _days(dt):
    """Get the days (floating point) from *d_t*.
    """
    return dt / np.timedelta64(1, 'D')


def gmst(utc_time):
    """Greenwich mean sidereal utc_time, in radians.

    As defined in the AIAA 2006 implementation:
    http://www.celestrak.com/publications/AIAA/2006-6753/
    """
    ut1 = jdays2000(utc_time) / 36525.0
    theta = 67310.54841 + ut1 * (876600 * 3600 + 8640184.812866 + ut1 *
                                 (0.093104 - ut1 * 6.2 * 10e-6))
    return np.deg2rad(theta / 240.0) % (2 * np.pi)


def _lmst(utc_time, longitude):
    """Local mean sidereal time, computed from *utc_time* and *longitude*.
    In radians.
    """
    return gmst(utc_time) + longitude


def sun_ecliptic_longitude(utc_time):
    """Ecliptic longitude of the sun at *utc_time*.
    """
    jdate = jdays2000(utc_time) / 36525.0
    # mean anomaly, rad
    m_a = np.deg2rad(357.52910 +
                     35999.05030 * jdate -
                     0.0001559 * jdate * jdate -
                     0.00000048 * jdate * jdate * jdate)
    # mean longitude, deg
    l_0 = 280.46645 + 36000.76983 * jdate + 0.0003032 * jdate * jdate
    d_l = ((1.914600 - 0.004817 * jdate - 0.000014 * jdate * jdate) * np.sin(m_a) +
           (0.019993 - 0.000101 * jdate) * np.sin(2 * m_a) + 0.000290 * np.sin(3 * m_a))
    # true longitude, deg
    l__ = l_0 + d_l
    return np.deg2rad(l__)


def sun_ra_dec(utc_time):
    """Right ascension and declination of the sun at *utc_time*.
    """
    jdate = jdays2000(utc_time) / 36525.0
    eps = np.deg2rad(23.0 + 26.0 / 60.0 + 21.448 / 3600.0 -
                     (46.8150 * jdate + 0.00059 * jdate * jdate -
                      0.001813 * jdate * jdate * jdate) / 3600)
    eclon = sun_ecliptic_longitude(utc_time)
    x__ = np.cos(eclon)
    y__ = np.cos(eps) * np.sin(eclon)
    z__ = np.sin(eps) * np.sin(eclon)
    r__ = np.sqrt(1.0 - z__ * z__)
    # sun declination
    declination = np.arctan2(z__, r__)
    # right ascension
    right_ascension = 2 * np.arctan2(y__, (x__ + r__))
    return right_ascension, declination


def _local_hour_angle(utc_time, longitude, right_ascension):
    """Hour angle at *utc_time* for the given *longitude* and
    *right_ascension*
    longitude in radians
    """
    return _lmst(utc_time, longitude) - right_ascension


def get_alt_az(utc_time, lon, lat):
    """Return sun altitude and azimuth from *utc_time*, *lon*, and *lat*.
    lon,lat in degrees
    The returned angles are given in radians.
    """
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    ra_, dec = sun_ra_dec(utc_time)
    h__ = _local_hour_angle(utc_time, lon, ra_)
    return (np.arcsin(np.sin(lat) * np.sin(dec) +
                      np.cos(lat) * np.cos(dec) * np.cos(h__)),
            np.arctan2(-np.sin(h__), (np.cos(lat) * np.tan(dec) -
                                      np.sin(lat) * np.cos(h__))))


def cos_zen(utc_time, lon, lat):
    """Cosine of the sun-zenith angle for *lon*, *lat* at *utc_time*.
    utc_time: datetime.datetime instance of the UTC time
    lon and lat in degrees.
    """
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    r_a, dec = sun_ra_dec(utc_time)
    h__ = _local_hour_angle(utc_time, lon, r_a)
    return (np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h__))


def sun_zenith_angle(utc_time, lon, lat):
    """Sun-zenith angle for *lon*, *lat* at *utc_time*.
    lon,lat in degrees.
    The angle returned is given in degrees
    """
    return np.rad2deg(np.arccos(cos_zen(utc_time, lon, lat)))


def sun_earth_distance_correction(utc_time):
    """Calculate the sun earth distance correction, relative to 1 AU.
    """
    year = 365.256363004
    # This is computed from
    # http://curious.astro.cornell.edu/question.php?number=582
    # AU = 149597870700.0
    # a = 149598261000.0
    # theta = (jdays2000(utc_time) - 2) * (2 * np.pi) / year
    # e = 0.01671123
    # r = a*(1-e*e)/(1+e * np.cos(theta))
    # corr_me = (r / AU) ** 2

    # from known software.
    corr = 1 - 0.0334 * np.cos(2 * np.pi * (jdays2000(utc_time) - 2) / year)

    return corr


def observer_position(time, lon, lat, alt):
    """Calculate observer ECI position.

    http://celestrak.com/columns/v02n03/
    """

    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    theta = (gmst(time) + lon) % (2 * np.pi)
    c = 1 / np.sqrt(1 + F * (F - 2) * np.sin(lat)**2)
    sq = c * (1 - F)**2

    achcp = (A * c + alt) * np.cos(lat)
    x = achcp * np.cos(theta)  # kilometers
    y = achcp * np.sin(theta)
    z = (A * sq + alt) * np.sin(lat)

    vx = -MFACTOR * y  # kilometers/second
    vy = MFACTOR * x
    vz = 0

    return (x, y, z), (vx, vy, vz)
