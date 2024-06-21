#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2011, 2013
#
# Author(s):
#
#   Martin Raspaud <martin.raspaud@smhi.se>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Angle and time-based astronomy functions.

Parts taken from http://www.geoastro.de/elevaz/basics/index.htm

Note on argument types
----------------------

Many of these functions accept Python datetime objects,
numpy datetime64 objects, or anything that can be turned
into a numpy array of datetime64 objects. These objects are inherently
64-bit so if other arguments (ex. longitude and latitude arrays) are
32-bit floats internal operations will be automatically promoted to
64-bit floating point numbers. Where possible these are then converted
back to 32-bit before being returned. In general scalar inputs will also
produce scalar outputs.

"""
from __future__ import annotations

import datetime
from typing import TypeAlias

import numpy as np
import numpy.typing as npt

from pyorbital import dt2np

F = 1 / 298.257223563  # Earth flattening WGS-84
A = 6378.137  # WGS84 Equatorial radius
MFACTOR = 7.292115E-5

ArrayOrFloat: TypeAlias = npt.ArrayLike | float
# numpy datetime or python datetime
DatetimeArrayLike: TypeAlias = npt.ArrayLike | np.datetime64 | datetime.datetime
TimedeltaArrayLike: TypeAlias = npt.ArrayLike | np.timedelta64


def jdays2000(utc_time: DatetimeArrayLike) -> np.ArrayLike[np.timedelta64]:
    """Get the days since year 2000.
    """
    return _days(dt2np(utc_time) - np.datetime64('2000-01-01T12:00'))


def jdays(utc_time: DatetimeArrayLike) -> float:
    """Get the julian day of *utc_time*.
    """
    return jdays2000(utc_time) + 2451545.0


def _days(dt: TimedeltaArrayLike) -> np.ArrayLike[np.float64]:
    """Get the days (floating point) from *d_t*.
    """
    if hasattr(dt, "shape"):
        dt = np.asanyarray(dt, dtype=np.timedelta64)
    return dt / np.timedelta64(1, 'D')


def gmst(utc_time: DatetimeArrayLike) -> npt.ArrayLike[np.float64]:
    """Greenwich mean sidereal utc_time, in radians.

    As defined in the AIAA 2006 implementation:
    http://www.celestrak.com/publications/AIAA/2006-6753/
    """
    ut1 = jdays2000(utc_time) / 36525.0
    theta = 67310.54841 + ut1 * (876600 * 3600 + 8640184.812866 + ut1 *
                                 (0.093104 - ut1 * 6.2 * 10e-6))
    return np.deg2rad(theta / 240.0) % (2 * np.pi)


def _lmst(utc_time: DatetimeArrayLike, longitude: ArrayOrFloat) -> npt.ArrayLike[np.float64]:
    """Local mean sidereal time, computed from *utc_time* and *longitude*.
    In radians.
    """
    return gmst(utc_time) + longitude


def sun_ecliptic_longitude(utc_time: datetime.datetime) -> npt.ArrayLike[np.float64]:
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


def sun_ra_dec(utc_time: datetime.datetime) -> tuple[npt.ArrayLike[np.float64], npt.ArrayLike[np.float64]]:
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


def _local_hour_angle(utc_time: DatetimeArrayLike, longitude: ArrayOrFloat, right_ascension: npt.ArrayLike[np.float64]) -> ArrayOrFloat:
    """Hour angle at *utc_time* for the given *longitude* and
    *right_ascension*
    longitude in radians
    """
    return _lmst(utc_time, longitude) - right_ascension


def get_alt_az(utc_time: DatetimeArrayLike, lon: ArrayOrFloat, lat: ArrayOrFloat) -> tuple[ArrayOrFloat, ArrayOrFloat]:
    """Return sun altitude and azimuth from *utc_time*, *lon*, and *lat*.

    lon,lat in degrees
    The returned angles are given in radians.
    """
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    ra_, dec = sun_ra_dec(utc_time)
    h__ = _local_hour_angle(utc_time, lon, ra_)
    alt_az = (np.arcsin(np.sin(lat) * np.sin(dec) +
                        np.cos(lat) * np.cos(dec) * np.cos(h__)),
              np.arctan2(-np.sin(h__), (np.cos(lat) * np.tan(dec) -
                                        np.sin(lat) * np.cos(h__))))
    if not isinstance(lon, float):
        alt_az = (alt_az[0].astype(lon.dtype), alt_az[1].astype(lon.dtype))
    return alt_az


def cos_zen(utc_time: DatetimeArrayLike, lon: ArrayOrFloat, lat: ArrayOrFloat) -> ArrayOrFloat:
    """Cosine of the sun-zenith angle for *lon*, *lat* at *utc_time*.
    utc_time: datetime.datetime instance of the UTC time
    lon and lat in degrees.
    """
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    r_a, dec = sun_ra_dec(utc_time)
    h__ = _local_hour_angle(utc_time, lon, r_a)
    csza = (np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h__))
    if not isinstance(lon, float):
        csza = csza.astype(lon.dtype)
    return csza


def sun_zenith_angle(utc_time: DatetimeArrayLike, lon: ArrayOrFloat, lat: ArrayOrFloat) -> ArrayOrFloat:
    """Sun-zenith angle for *lon*, *lat* at *utc_time*.
    lon,lat in degrees.
    The angle returned is given in degrees
    """
    sza = np.rad2deg(np.arccos(cos_zen(utc_time, lon, lat)))
    if not isinstance(lon, float):
        sza = sza.astype(lon.dtype)
    return sza


def sun_earth_distance_correction(utc_time: DatetimeArrayLike) -> ArrayOrFloat:
    """Calculate the sun earth distance correction, relative to 1 AU.
    """
    # Computation according to
    #  https://web.archive.org/web/20150117190838/http://curious.astro.cornell.edu/question.php?number=582
    # with
    #  Astronomical unit: AU = 149597870700.0 meter (https://ssd.jpl.nasa.gov/glossary/au.html)
    #  Semi-major axis:    a = 1.00000261 AU = 149598261150.0 m (https://ssd.jpl.nasa.gov/planets/approx_pos.html)
    #  Eccentricity:       e = 0.01671123 (https://ssd.jpl.nasa.gov/planets/approx_pos.html)
    #  Length of year:  year = 365.25636 days (https://ssd.jpl.nasa.gov/astro_par.html)
    #  Perihelion:         p = 3 (day of year with Earth in perihelion, varies between Jan 2 and Jan 5,
    #                             http://www.astropixels.com/ephemeris/perap2001.html
    #                             https://web.archive.org/web/20080328044924/http://aa.usno.navy.mil/data/docs/EarthSeasons)
    # Formula:
    #  theta = (jdays2000(utc_time) - p) * (2 * np.pi) / year
    #  r = a * (1 - e * e) / (1 + e * np.cos(theta))
    #  corr := r/AU
    #        =  a * (1 - e * e) / AU / (1 + e * np.cos(theta))
    #       "=" a * (1 - e * e) / AU * (1 - e * np.cos(theta))
    #       "=" 1 - 0.0167 * np.cos(theta)

    corr = 1 - 0.0167 * np.cos(2 * np.pi * (jdays2000(utc_time) - 3) / 365.25636)
    return corr


def observer_position(
        utc_time: DatetimeArrayLike, lon: ArrayOrFloat, lat: ArrayOrFloat, alt: ArrayOrFloat
) -> tuple[tuple[ArrayOrFloat, ArrayOrFloat, ArrayOrFloat], tuple[ArrayOrFloat, ArrayOrFloat, ArrayOrFloat]]:
    """Calculate observer ECI position.

    http://celestrak.com/columns/v02n03/
    """

    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    theta = (gmst(utc_time) + lon) % (2 * np.pi)
    c = 1 / np.sqrt(1 + F * (F - 2) * np.sin(lat)**2)
    sq = c * (1 - F)**2

    achcp = (A * c + alt) * np.cos(lat)
    x = achcp * np.cos(theta)  # kilometers
    y = achcp * np.sin(theta)
    z = (A * sq + alt) * np.sin(lat)

    vx = -MFACTOR * y  # kilometers/second
    vy = MFACTOR * x
    vz = _float_to_sibling_result(0.0, vx)

    if not isinstance(lon, float):
        x = x.astype(lon.dtype, copy=False)
        y = y.astype(lon.dtype, copy=False)
        z = z.astype(lon.dtype, copy=False)
        vx = vx.astype(lon.dtype, copy=False)
        vy = vy.astype(lon.dtype, copy=False)
        vz = vz.astype(lon.dtype, copy=False)  # type: ignore[union-attr]
    return (x, y, z), (vx, vy, vz)


def _float_to_sibling_result(
        result_to_convert: float,
        template_result: ArrayOrFloat,
) -> ArrayOrFloat:
    """Convert a scalar to the same type as another return type.

    This is mostly used to make a static value consistent with the types of
    other returned values.

    """
    if isinstance(template_result, float):
        return result_to_convert
    # get any array like object that might be wrapped by our template (ex. xarray DataArray)
    array_like = template_result if hasattr(template_result, "__array_function__") else template_result.data
    array_convert = np.asarray(result_to_convert, like=array_like)
    if not hasattr(template_result, "__array_function__"):
        # the template result has some wrapper class (likely xarray DataArray)
        # recreate the wrapper object
        array_convert = template_result.__class__(array_convert)
    return array_convert
