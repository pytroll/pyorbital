#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012, 2013, 2014, 2015.

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

"""Module to compute geolocalization of a satellite scene.
"""

# TODO:
# - Attitude correction
# - project on an ellipsoid instead of a sphere
# - optimize !!!
# - test !!!

from __future__ import print_function

import numpy as np
from numpy import cos, sin, sqrt
from datetime import timedelta
from pyorbital.orbital import Orbital

a = 6378.137  # km
b = 6356.75231414  # km, GRS80
# b = 6356.752314245 # km, WGS84


def geodetic_lat(point, a=a, b=b):

    x, y, z = point
    r = np.sqrt(x * x + y * y)
    geoc_lat = np.arctan2(z, r)

    geod_lat = geoc_lat
    e2 = (a * a - b * b) / (a * a)
    while True:
        phi = geod_lat
        C = 1 / sqrt(1 - e2 * sin(phi)**2)
        geod_lat = np.arctan2(z + a * C * e2 * sin(phi), r)
        if np.allclose(geod_lat, phi):
            return geod_lat


def subpoint(query_point, a=a, b=b):
    """Get the point on the ellipsoid under the *query_point*.
    """
    x, y, z = query_point
    r = sqrt(x * x + y * y)

    lat = geodetic_lat(query_point)
    lon = np.arctan2(y, x)
    e2_ = (a * a - b * b) / (a * a)
    n__ = a / sqrt(1 - e2_ * sin(lat)**2)
    nx_ = n__ * cos(lat) * cos(lon)
    ny_ = n__ * cos(lat) * sin(lon)
    nz_ = (1 - e2_) * n__ * sin(lat)

    return np.vstack([nx_, ny_, nz_])


class ScanGeometry(object):

    """Description of the geometry of an instrument.

    *fovs* is the x and y viewing angles of the instrument. y is zero if the we
    talk about scanlines of course. *times* is the time of viewing of each
    angle relative to the start of the scanning, so it should have the same
    size as the *fovs*. *attitude* is the attitude correction to apply.
    """

    def __init__(self,
                 fovs,
                 times,
                 attitude=(0, 0, 0)):
        self.fovs = np.array(fovs)
        self._times = np.array(times)
        self.attitude = attitude

    def vectors(self, pos, vel, roll=0.0, pitch=0.0, yaw=0.0):
        """Get unit vectors pointing to the different pixels.

        *pos* and *vel* are column vectors, or matrices of column
        vectors. Returns vectors as stacked rows.
        """
        # TODO: yaw steering mode !

        # Fake nadir: This is the intersection point between the satellite
        # looking down at the centre of the ellipsoid and the surface of the
        # ellipsoid. Nadir on the other hand is the point which vertical goes
        # through the satellite...
        #nadir = -pos / vnorm(pos)

        nadir = subpoint(-pos)
        nadir /= vnorm(nadir)

        # x is along track (roll)
        x = vel / vnorm(vel)

        # y is cross track (pitch)
        y = np.cross(nadir, vel, 0, 0, 0)
        y /= vnorm(y)

        # rotate first around x
        x_rotated = qrotate(nadir, x, self.fovs[:, 0] + roll)
        # then around y
        xy_rotated = qrotate(x_rotated, y, self.fovs[:, 1] + pitch)
        # then around z
        return qrotate(xy_rotated, nadir, yaw)

    def times(self, start_of_scan):
        tds = [timedelta(seconds=i) for i in self._times]
        return np.array(tds) + start_of_scan


class Quaternion(object):

    def __init__(self, scalar, vector):
        self.__x, self.__y, self.__z = vector
        self.__w = scalar

    def rotation_matrix(self):
        x, y, z, w = self.__x, self.__y, self.__z, self.__w
        zero = np.zeros_like(x)
        return np.array(
            ((w**2 + x**2 - y**2 - z**2,
              2 * x * y + 2 * z * w,
              2 * x * z - 2 * y * w,
              zero),
             (2 * x * y - 2 * z * w,
              w**2 - x**2 + y**2 - z**2,
              2 * y * z + 2 * x * w,
              zero),
             (2 * x * z + 2 * y * w,
              2 * y * z - 2 * x * w,
              w**2 - x**2 - y**2 + z**2,
              zero),
             (zero, zero, zero, w**2 + x**2 + y**2 + z**2)))


def qrotate(vector, axis, angle):
    """Rotate *vector* around *axis* by *angle* (in radians).

    *vector* is a matrix of column vectors, as is *axis*.
    This function uses quaternion rotation.
    """
    n_axis = axis / vnorm(axis)
    sin_angle = np.expand_dims(sin(angle / 2), 0)
    if np.rank(n_axis) == 1:
        n_axis = np.expand_dims(n_axis, 1)
        p__ = np.dot(n_axis, sin_angle)[:, np.newaxis]
    else:
        p__ = n_axis * sin_angle

    q__ = Quaternion(cos(angle / 2), p__)
    return np.einsum("kj, ikj->ij",
                     vector,
                     q__.rotation_matrix()[:3, :3])


# DIRTY STUFF. Needed the get_lonlatalt function to work on pos directly if
# we want to print out lonlats in the end.
from pyorbital import astronomy
from pyorbital.orbital import *


def get_lonlatalt(pos, utc_time):
    """Calculate sublon, sublat and altitude of satellite, considering the
    earth an ellipsoid.

    http://celestrak.com/columns/v02n03/
    """
    (pos_x, pos_y, pos_z) = pos / XKMPER
    lon = ((np.arctan2(pos_y * XKMPER, pos_x * XKMPER) - astronomy.gmst(utc_time))
           % (2 * np.pi))

    lon = np.where(lon > np.pi, lon - np.pi * 2, lon)
    lon = np.where(lon <= -np.pi, lon + np.pi * 2, lon)

    r = np.sqrt(pos_x ** 2 + pos_y ** 2)
    lat = np.arctan2(pos_z, r)
    e2 = F * (2 - F)

    while True:
        lat2 = lat
        c = 1 / (np.sqrt(1 - e2 * (np.sin(lat2) ** 2)))
        lat = np.arctan2(pos_z + c * e2 * np.sin(lat2), r)
        if np.all(abs(lat - lat2) < 1e-10):
            break
    alt = r / np.cos(lat) - c
    alt *= A
    return np.rad2deg(lon), np.rad2deg(lat), alt

# END OF DIRTY STUFF


def compute_pixels(tle, sgeom, times, rpy=(0.0, 0.0, 0.0)):
    """Compute cartesian coordinates of the pixels in instrument scan.
    """
    (tle1, tle2) = tle
    orb = Orbital("mysatellite", line1=tle1, line2=tle2)

    # get position and velocity for each time of each pixel
    pos, vel = orb.get_position(times, normalize=False)

    # now, get the vectors pointing to each pixel
    vectors = sgeom.vectors(pos, vel, *rpy)

    # compute intersection of lines (directed by vectors and passing through
    # (0, 0, 0)) and ellipsoid. Derived from:
    # http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    # do the computation between line and ellipsoid (WGS 84)
    # NB: AAPP uses GRS 80...
    centre = -pos
    a__ = 6378.137  # km
    # b__ = 6356.75231414 # km, GRS80
    b__ = 6356.752314245  # km, WGS84
    radius = np.array([[1 / a__, 1 / a__, 1 / b__]]).T
    xr_ = vectors * radius
    cr_ = centre * radius
    ldotc = np.einsum("ij,ij->j", xr_, cr_)
    lsq = np.einsum("ij,ij->j", xr_, xr_)
    csq = np.einsum("ij,ij->j", cr_, cr_)

    d1_ = (ldotc - np.sqrt(ldotc ** 2 - csq * lsq + lsq)) / lsq

    # return the actual pixel positions
    return vectors * d1_ - centre


def norm(v):
    return np.sqrt(np.dot(v, v.conj()))


def mnorm(m, axis=None):
    """norm of a matrix of vectors stacked along the *axis* dimension.
    """
    if axis is None:
        axis = np.rank(m) - 1
    return np.sqrt((m**2).sum(axis))


def vnorm(m):
    """norms of a matrix of column vectors.
    """
    return np.sqrt((m**2).sum(0))


def hnorm(m):
    """norms of a matrix of row vectors.
    """
    return np.sqrt((m**2).sum(1))

if __name__ == '__main__':
    # NOAA 18 (from the 2011-10-12, 16:55 utc)
    # 1 28654U 05018A   11284.35271227  .00000478  00000-0  28778-3 0  9246
    # 2 28654  99.0096 235.8581 0014859 135.4286 224.8087 14.11526826329313

    noaa18_tle1 = "1 28654U 05018A   11284.35271227  .00000478  00000-0  28778-3 0  9246"
    noaa18_tle2 = "2 28654  99.0096 235.8581 0014859 135.4286 224.8087 14.11526826329313"

    from datetime import datetime
    t = datetime(2011, 10, 12, 13, 45)

    # edge and centre of an avhrr scanline
    # sgeom = ScanGeometry([(-0.9664123687741623, 0),
    #                      (0, 0)],
    #                     [0, 0.0, ])
    # print compute_pixels((noaa18_tle1, noaa18_tle2), sgeom, t)

    # avhrr swath
    scanline_nb = 1

    # building the avhrr angles, 2048 pixels from +55.37 to -55.37 degrees
    avhrr = np.vstack(((np.arange(2048) - 1023.5) / 1024 * np.deg2rad(-55.37),
                       np.zeros((2048,)))).transpose()
    avhrr = np.tile(avhrr, [scanline_nb, 1])
    # building the corresponding times array
    offset = np.arange(scanline_nb) * 0.1667
    times = (np.tile(np.arange(2048) * 0.000025 + 0.0025415, [scanline_nb, 1])
             + np.expand_dims(offset, 1))
    # build the scan geometry object
    sgeom = ScanGeometry(avhrr, times.ravel())

    # print the lonlats for the pixel positions
    s_times = sgeom.times(t)
    pixels_pos = compute_pixels((noaa18_tle1, noaa18_tle2), sgeom, s_times)
    print(get_lonlatalt(pixels_pos, s_times))
