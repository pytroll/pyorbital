#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012, 2013, 2014, 2015.

# Author(s):

#   Esben S. Nielsen <esn@dmi.dk>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
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

"""Module for computing the orbital parameters of satellites."""

import logging
import warnings
from datetime import datetime, timedelta

import numpy as np
from scipy import optimize

from pyorbital import astronomy, dt2np, tlefile

try:
    import dask.array as da
    has_dask = True
except ImportError:
    da = None
    has_dask = False

try:
    import xarray as xr
    has_xarray = True
except ImportError:
    xr = None
    has_xarray = False

logger = logging.getLogger(__name__)

ECC_EPS = 1.0e-6  # Too low for computing further drops.
ECC_LIMIT_LOW = -1.0e-3
ECC_LIMIT_HIGH = 1.0 - ECC_EPS  # Too close to 1
ECC_ALL = 1.0e-4

EPS_COS = 1.5e-12

NR_EPS = 1.0e-12

CK2 = 5.413080e-4
CK4 = 0.62098875e-6
E6A = 1.0e-6
QOMS2T = 1.88027916e-9
S = 1.01222928
S0 = 78.0
XJ3 = -0.253881e-5
XKE = 0.743669161e-1
XKMPER = 6378.135
XMNPDA = 1440.0
# MFACTOR = 7.292115E-5
AE = 1.0
SECDAY = 8.6400E4

F = 1 / 298.257223563  # Earth flattening WGS-84
A = 6378.137  # WGS84 Equatorial radius


SGDP4_ZERO_ECC = 0
SGDP4_DEEP_NORM = 1
SGDP4_NEAR_SIMP = 2
SGDP4_NEAR_NORM = 3

KS = AE * (1.0 + S0 / XKMPER)
A3OVK2 = (-XJ3 / CK2) * AE**3


class OrbitalError(Exception):
    pass


def get_observer_look(sat_lon, sat_lat, sat_alt, utc_time, lon, lat, alt):
    """Calculate observers look angle to a satellite.
    http://celestrak.com/columns/v02n02/

    :utc_time: Observation time (datetime object)
    :lon: Longitude of observer position on ground in degrees east
    :lat: Latitude of observer position on ground in degrees north
    :alt: Altitude above sea-level (geoid) of observer position on ground in km

    :return: (Azimuth, Elevation)
    """
    (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = astronomy.observer_position(
        utc_time, sat_lon, sat_lat, sat_alt)

    (opos_x, opos_y, opos_z), (ovel_x, ovel_y, ovel_z) = \
        astronomy.observer_position(utc_time, lon, lat, alt)

    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    theta = (astronomy.gmst(utc_time) + lon) % (2 * np.pi)

    rx = pos_x - opos_x
    ry = pos_y - opos_y
    rz = pos_z - opos_z

    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    top_s = sin_lat * cos_theta * rx + \
        sin_lat * sin_theta * ry - cos_lat * rz
    top_e = -sin_theta * rx + cos_theta * ry
    top_z = cos_lat * cos_theta * rx + \
        cos_lat * sin_theta * ry + sin_lat * rz

    # Azimuth is undefined when elevation is 90 degrees, 180 (pi) will be returned.
    az_ = np.arctan2(-top_e, top_s) + np.pi
    az_ = np.mod(az_, 2 * np.pi)  # Needed on some platforms

    rg_ = np.sqrt(rx * rx + ry * ry + rz * rz)

    top_z_divided_by_rg_ = top_z / rg_

    # Due to rounding top_z can be larger than rg_ (when el_ ~ 90).
    top_z_divided_by_rg_ = top_z_divided_by_rg_.clip(max=1)
    el_ = np.arcsin(top_z_divided_by_rg_)

    return np.rad2deg(az_), np.rad2deg(el_)


class Orbital(object):

    """Class for orbital computations.

    The *satellite* parameter is the name of the satellite to work on and is
    used to retrieve the right TLE data for internet or from *tle_file* in case
    it is provided.
    """

    def __init__(self, satellite, tle_file=None, line1=None, line2=None):
        satellite = satellite.upper()
        self.satellite_name = satellite
        self.tle = tlefile.read(satellite, tle_file=tle_file,
                                line1=line1, line2=line2)
        self.orbit_elements = OrbitElements(self.tle)
        self._sgdp4 = _SGDP4(self.orbit_elements)

    def __str__(self):
        return self.satellite_name + " " + str(self.tle)

    def get_last_an_time(self, utc_time):
        """Calculate time of last ascending node relative to the
        specified time
        """

        # Propagate backwards to ascending node
        dt = np.timedelta64(10, 'm')
        t_old = utc_time
        t_new = t_old - dt
        pos0, vel0 = self.get_position(t_old, normalize=False)
        pos1, vel1 = self.get_position(t_new, normalize=False)
        while not (pos0[2] > 0 and pos1[2] < 0):
            pos0 = pos1
            t_old = t_new
            t_new = t_old - dt
            pos1, vel1 = self.get_position(t_new, normalize=False)

        # Return if z within 1 km of an
        if np.abs(pos0[2]) < 1:
            return t_old
        elif np.abs(pos1[2]) < 1:
            return t_new

        # Bisect to z within 1 km
        while np.abs(pos1[2]) > 1:
            # pos0, vel0 = pos1, vel1
            dt = (t_old - t_new) / 2
            t_mid = t_old - dt
            pos1, vel1 = self.get_position(t_mid, normalize=False)
            if pos1[2] > 0:
                t_old = t_mid
            else:
                t_new = t_mid

        return t_mid

    def get_position(self, utc_time, normalize=True):
        """Get the cartesian position and velocity from the satellite."""
        kep = self._sgdp4.propagate(utc_time)
        pos, vel = kep2xyz(kep)

        if normalize:
            pos /= XKMPER
            vel /= XKMPER * XMNPDA / SECDAY

        return pos, vel

    def get_lonlatalt(self, utc_time):
        """Calculate sublon, sublat and altitude of satellite.

        http://celestrak.com/columns/v02n03/
        """
        (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = self.get_position(
            utc_time, normalize=True)

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

    def find_aos(self, utc_time, lon, lat):
        pass

    def find_aol(self, utc_time, lon, lat):
        pass

    def get_observer_look(self, utc_time, lon, lat, alt):
        """Calculate observers look angle to a satellite.
        http://celestrak.com/columns/v02n02/

        utc_time: Observation time (datetime object)
        lon: Longitude of observer position on ground in degrees east
        lat: Latitude of observer position on ground in degrees north
        alt: Altitude above sea-level (geoid) of observer position on ground in km

        Return: (Azimuth, Elevation)
        """

        utc_time = dt2np(utc_time)
        (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = self.get_position(
            utc_time, normalize=False)
        (opos_x, opos_y, opos_z), (ovel_x, ovel_y, ovel_z) = \
            astronomy.observer_position(utc_time, lon, lat, alt)

        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)

        theta = (astronomy.gmst(utc_time) + lon) % (2 * np.pi)

        rx = pos_x - opos_x
        ry = pos_y - opos_y
        rz = pos_z - opos_z

        sin_lat = np.sin(lat)
        cos_lat = np.cos(lat)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        top_s = sin_lat * cos_theta * rx + \
            sin_lat * sin_theta * ry - cos_lat * rz
        top_e = -sin_theta * rx + cos_theta * ry
        top_z = cos_lat * cos_theta * rx + \
            cos_lat * sin_theta * ry + sin_lat * rz

        az_ = np.arctan(-top_e / top_s)

        az_ = np.where(top_s > 0, az_ + np.pi, az_)
        az_ = np.where(az_ < 0, az_ + 2 * np.pi, az_)

        rg_ = np.sqrt(rx * rx + ry * ry + rz * rz)
        el_ = np.arcsin(top_z / rg_)

        return np.rad2deg(az_), np.rad2deg(el_)

    def get_orbit_number(self, utc_time, tbus_style=False, as_float=False):
        """Calculate orbit number at specified time.

        Args:
            tbus_style: If True, use TBUS-style orbit numbering (TLE orbit number + 1)
            as_float: Return a continuous orbit number as float.
        """
        utc_time = np.datetime64(utc_time)
        try:
            dt = astronomy._days(utc_time - self.orbit_elements.an_time)
            orbit_period = astronomy._days(self.orbit_elements.an_period)
        except AttributeError:
            pos_epoch, vel_epoch = self.get_position(self.tle.epoch,
                                                     normalize=False)
            if np.abs(pos_epoch[2]) > 1 or not vel_epoch[2] > 0:
                # Epoch not at ascending node
                self.orbit_elements.an_time = self.get_last_an_time(
                    self.tle.epoch)
            else:
                # Epoch at ascending node (z < 1 km) and positive v_z
                self.orbit_elements.an_time = self.tle.epoch

            self.orbit_elements.an_period = self.orbit_elements.an_time - \
                self.get_last_an_time(self.orbit_elements.an_time
                                      - np.timedelta64(10, 'm'))

            dt = astronomy._days(utc_time - self.orbit_elements.an_time)
            orbit_period = astronomy._days(self.orbit_elements.an_period)

        orbit = self.tle.orbit + dt / orbit_period + \
            self.tle.mean_motion_derivative * dt ** 2 + \
            self.tle.mean_motion_sec_derivative * dt ** 3
        if not as_float:
            orbit = int(orbit)

        if tbus_style:
            orbit += 1

        return orbit

    def get_next_passes(self, utc_time, length, lon, lat, alt, tol=0.001, horizon=0):
        """Calculate passes for the next hours for a given start time and a
        given observer.

        Original by Martin.

        :utc_time: Observation time (datetime object)
        :length: Number of hours to find passes (int)
        :lon: Longitude of observer position on ground (float)
        :lat: Latitude of observer position on ground (float)
        :alt: Altitude above sea-level (geoid) of observer position on ground (float)
        :tol: precision of the result in seconds
        :horizon: the elevation of horizon to compute risetime and falltime.

        :return: [(rise-time, fall-time, max-elevation-time), ...]
        """

        def elevation(minutes):
            """Compute the elevation."""
            return self.get_observer_look(utc_time +
                                          timedelta(
                                              minutes=np.float64(minutes)),
                                          lon, lat, alt)[1] - horizon

        def elevation_inv(minutes):
            """Compute the inverse of elevation."""
            return -elevation(minutes)

        def get_root(fun, start, end, tol=0.01):
            """Root finding scheme"""
            x_0 = end
            x_1 = start
            fx_0 = fun(end)
            fx_1 = fun(start)
            if abs(fx_0) < abs(fx_1):
                fx_0, fx_1 = fx_1, fx_0
                x_0, x_1 = x_1, x_0

            x_n = optimize.brentq(fun, x_0, x_1)
            return x_n

        def get_max_parab(fun, start, end, tol=0.01):
            """Successive parabolic interpolation."""
            a = float(start)
            c = float(end)
            b = (a + c) / 2.0

            f_a = fun(a)
            f_b = fun(b)
            f_c = fun(c)

            x = b
            with np.errstate(invalid='raise'):
                while True:
                    try:
                        x = x - 0.5 * (((b - a) ** 2 * (f_b - f_c)
                                        - (b - c) ** 2 * (f_b - f_a)) /
                                       ((b - a) * (f_b - f_c) - (b - c) * (f_b - f_a)))
                    except FloatingPointError:
                        return b
                    if abs(b - x) <= tol:
                        return x
                    f_x = fun(x)
                    # sometimes the estimation diverges... return best guess
                    if f_x > f_b:
                        logger.info('Parabolic interpolation did not converge, returning best guess so far.')
                        return b

                    a, b, c = (a + x) / 2.0, x, (x + c) / 2.0
                    f_a, f_b, f_c = fun(a), f_x, fun(c)

        # every minute
        times = utc_time + np.array([timedelta(minutes=minutes)
                                     for minutes in range(length * 60)])
        elev = self.get_observer_look(times, lon, lat, alt)[1] - horizon
        zcs = np.where(np.diff(np.sign(elev)))[0]
        res = []
        risetime = None
        for guess in zcs:
            horizon_mins = get_root(
                elevation, guess, guess + 1.0, tol=tol / 60.0)
            horizon_time = utc_time + timedelta(minutes=horizon_mins)
            if elev[guess] < 0:
                risetime = horizon_time
                risemins = horizon_mins
            else:
                falltime = horizon_time
                fallmins = horizon_mins
                if risetime:
                    int_start = max(0, int(np.floor(risemins)))
                    int_end = min(len(elev), int(np.ceil(fallmins) + 1))
                    middle = int_start + np.argmax(elev[int_start:int_end])
                    highest = utc_time + \
                        timedelta(minutes=get_max_parab(
                            elevation_inv,
                            max(risemins, middle - 1), min(fallmins, middle + 1),
                            tol=tol / 60.0
                        ))
                    res += [(risetime, falltime, highest)]
                risetime = None
        return res

    def _get_time_at_horizon(self, utc_time, obslon, obslat, **kwargs):
        """Get the time closest in time to *utc_time* when the
        satellite is at the horizon relative to the position of an observer on
        ground (altitude = 0)

        Note: This is considered deprecated and it's functionality is currently
        replaced by 'get_next_passes'.
        """
        warnings.warn("_get_time_at_horizon is replaced with get_next_passes",
                      DeprecationWarning)
        if "precision" in kwargs:
            precision = kwargs['precision']
        else:
            precision = timedelta(seconds=0.001)
        if "max_iterations" in kwargs:
            nmax_iter = kwargs["max_iterations"]
        else:
            nmax_iter = 100

        sec_step = 0.5
        t_step = timedelta(seconds=sec_step / 2.0)

        # Local derivative:
        def fprime(timex):
            el0 = self.get_observer_look(timex - t_step,
                                         obslon, obslat, 0.0)[1]
            el1 = self.get_observer_look(timex + t_step,
                                         obslon, obslat, 0.0)[1]
            return el0, (abs(el1) - abs(el0)) / sec_step

        tx0 = utc_time - timedelta(seconds=1.0)
        tx1 = utc_time
        idx = 0
        # eps = 500.
        eps = 100.
        while abs(tx1 - tx0) > precision and idx < nmax_iter:
            tx0 = tx1
            fpr = fprime(tx0)
            # When the elevation is high the scale is high, and when
            # the elevation is low the scale is low
            # var_scale = np.abs(np.sin(fpr[0] * np.pi/180.))
            # var_scale = np.sqrt(var_scale)
            var_scale = np.abs(fpr[0])
            tx1 = tx0 - timedelta(seconds=(eps * var_scale * fpr[1]))
            idx = idx + 1
            # print idx, tx0, tx1, var_scale, fpr
            if abs(tx1 - utc_time) < precision and idx < 2:
                tx1 = tx1 + timedelta(seconds=1.0)

        if abs(tx1 - tx0) <= precision and idx < nmax_iter:
            return tx1
        else:
            return None

    def utc2local(self, utc_time):
        """Convert UTC to local time."""
        lon, _, _ = self.get_lonlatalt(utc_time)
        return utc_time + timedelta(hours=lon * 24 / 360.0)

    def get_equatorial_crossing_time(self, tstart, tend, node='ascending', local_time=False,
                                     rtol=1E-9):
        """Estimate the equatorial crossing time of an orbit.

        The crossing time is determined via the orbit number, which increases by one if the
        spacecraft passes the ascending node at the equator. A bisection algorithm is used to find
        the time of that passage.

        Args:
            tstart: Start time of the orbit
            tend: End time of the orbit. Orbit number at the end must be at least one greater than
                at the start. If there are multiple revolutions in the given time interval, the
                crossing time of the last revolution in that interval will be computed.
            node: Specifies whether to compute the crossing time at the ascending or descending
                node. Choices: ('ascending', 'descending').
            local_time: By default the UTC crossing time is returned. Use this flag to convert UTC
                to local time.
            rtol: Tolerance of the bisection algorithm. The smaller the tolerance, the more accurate
                the result.
        """
        # Determine orbit number at the start and end of the orbit.
        n_start = self.get_orbit_number(tstart, as_float=True)
        n_end = self.get_orbit_number(tend, as_float=True)
        if int(n_end) - int(n_start) == 0:
            # Orbit doesn't cross the equator in the given time interval
            return None
        elif n_end - n_start > 1:
            warnings.warn('Multiple revolutions between start and end time. Computing crossing '
                          'time for the last revolution in that interval.')

        # Let n'(t) = n(t) - offset. Determine offset so that n'(tstart) < 0 and n'(tend) > 0 and
        # n'(tcross) = 0.
        offset = int(n_end)
        if node == 'descending':
            offset = offset + 0.5

        # Use bisection algorithm to find the root of n'(t), which is the crossing time. The
        # algorithm requires continuous time coordinates, so convert timestamps to microseconds
        # since 1970.
        time_unit = 'us'  # same precision as datetime

        def _nprime(time_f):
            """Continuous orbit number as a function of time."""
            time64 = np.datetime64(int(time_f), time_unit)
            n = self.get_orbit_number(time64, as_float=True)
            return n - offset

        try:
            tcross = optimize.bisect(_nprime,
                                     a=np.datetime64(tstart, time_unit).astype(np.int64),
                                     b=np.datetime64(tend, time_unit).astype(np.int64),
                                     rtol=rtol)
        except ValueError:
            # Bisection did not converge
            return None
        tcross = np.datetime64(int(tcross), time_unit).astype(datetime)

        # Convert UTC to local time
        if local_time:
            tcross = self.utc2local(tcross)

        return tcross


class OrbitElements(object):

    """Class holding the orbital elements.
    """

    def __init__(self, tle):
        self.epoch = tle.epoch
        self.excentricity = tle.excentricity
        self.inclination = np.deg2rad(tle.inclination)
        self.right_ascension = np.deg2rad(tle.right_ascension)
        self.arg_perigee = np.deg2rad(tle.arg_perigee)
        self.mean_anomaly = np.deg2rad(tle.mean_anomaly)

        self.mean_motion = tle.mean_motion * (np.pi * 2 / XMNPDA)
        self.mean_motion_derivative = tle.mean_motion_derivative * \
            np.pi * 2 / XMNPDA ** 2
        self.mean_motion_sec_derivative = tle.mean_motion_sec_derivative * \
            np.pi * 2 / XMNPDA ** 3
        self.bstar = tle.bstar * AE

        n_0 = self.mean_motion
        k_e = XKE
        k_2 = CK2
        i_0 = self.inclination
        e_0 = self.excentricity

        a_1 = (k_e / n_0) ** (2.0 / 3)
        delta_1 = ((3 / 2.0) * (k_2 / a_1**2) * ((3 * np.cos(i_0)**2 - 1) /
                                                 (1 - e_0**2)**(2.0 / 3)))

        a_0 = a_1 * (1 - delta_1 / 3 - delta_1**2 - (134.0 / 81) * delta_1**3)

        delta_0 = ((3 / 2.0) * (k_2 / a_0**2) * ((3 * np.cos(i_0)**2 - 1) /
                                                 (1 - e_0**2)**(2.0 / 3)))

        # original mean motion
        n_0pp = n_0 / (1 + delta_0)
        self.original_mean_motion = n_0pp

        # semi major axis
        a_0pp = a_0 / (1 - delta_0)
        self.semi_major_axis = a_0pp

        self.period = np.pi * 2 / n_0pp

        self.perigee = (a_0pp * (1 - e_0) / AE - AE) * XKMPER

        self.right_ascension_lon = (self.right_ascension
                                    - astronomy.gmst(self.epoch))

        if self.right_ascension_lon > np.pi:
            self.right_ascension_lon -= 2 * np.pi


class _SGDP4(object):

    """Class for the SGDP4 computations.
    """

    def __init__(self, orbit_elements):
        self.mode = None

        # perigee = orbit_elements.perigee
        self.eo = orbit_elements.excentricity
        self.xincl = orbit_elements.inclination
        self.xno = orbit_elements.original_mean_motion
        # k_2 = CK2
        # k_4 = CK4
        # k_e = XKE
        self.bstar = orbit_elements.bstar
        self.omegao = orbit_elements.arg_perigee
        self.xmo = orbit_elements.mean_anomaly
        self.xnodeo = orbit_elements.right_ascension
        self.t_0 = orbit_elements.epoch
        self.xn_0 = orbit_elements.mean_motion
        # A30 = -XJ3 * AE**3

        if not(0 < self.eo < ECC_LIMIT_HIGH):
            raise OrbitalError('Eccentricity out of range: %e' % self.eo)
        elif not((0.0035 * 2 * np.pi / XMNPDA) < self.xn_0 < (18 * 2 * np.pi / XMNPDA)):
            raise OrbitalError('Mean motion out of range: %e' % self.xn_0)
        elif not(0 < self.xincl < np.pi):
            raise OrbitalError('Inclination out of range: %e' % self.xincl)

        if self.eo < 0:
            self.mode = self.SGDP4_ZERO_ECC
            return

        self.cosIO = np.cos(self.xincl)
        self.sinIO = np.sin(self.xincl)
        theta2 = self.cosIO**2
        theta4 = theta2 ** 2
        self.x3thm1 = 3.0 * theta2 - 1.0
        self.x1mth2 = 1.0 - theta2
        self.x7thm1 = 7.0 * theta2 - 1.0

        a1 = (XKE / self.xn_0) ** (2. / 3)
        betao2 = 1.0 - self.eo**2
        betao = np.sqrt(betao2)
        temp0 = 1.5 * CK2 * self.x3thm1 / (betao * betao2)
        del1 = temp0 / (a1**2)
        a0 = a1 * \
            (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)))
        del0 = temp0 / (a0**2)
        self.xnodp = self.xn_0 / (1.0 + del0)
        self.aodp = (a0 / (1.0 - del0))
        self.perigee = (self.aodp * (1.0 - self.eo) - AE) * XKMPER
        self.apogee = (self.aodp * (1.0 + self.eo) - AE) * XKMPER
        self.period = (2 * np.pi * 1440.0 / XMNPDA) / self.xnodp

        if self.period >= 225:
            # Deep-Space model
            self.mode = SGDP4_DEEP_NORM
        elif self.perigee < 220:
            # Near-space, simplified equations
            self.mode = SGDP4_NEAR_SIMP
        else:
            # Near-space, normal equations
            self.mode = SGDP4_NEAR_NORM

        if self.perigee < 156:
            s4 = self.perigee - 78
            if s4 < 20:
                s4 = 20

            qoms24 = ((120 - s4) * (AE / XKMPER))**4
            s4 = (s4 / XKMPER + AE)
        else:
            s4 = KS
            qoms24 = QOMS2T

        pinvsq = 1.0 / (self.aodp**2 * betao2**2)
        tsi = 1.0 / (self.aodp - s4)
        self.eta = self.aodp * self.eo * tsi
        etasq = self.eta**2
        eeta = self.eo * self.eta
        psisq = np.abs(1.0 - etasq)
        coef = qoms24 * tsi**4
        coef_1 = coef / psisq**3.5

        self.c2 = (coef_1 * self.xnodp * (self.aodp *
                                          (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                                          (0.75 * CK2) * tsi / psisq * self.x3thm1 *
                                          (8.0 + 3.0 * etasq * (8.0 + etasq))))

        self.c1 = self.bstar * self.c2

        self.c4 = (2.0 * self.xnodp * coef_1 * self.aodp * betao2 * (
            self.eta * (2.0 + 0.5 * etasq) + self.eo * (0.5 + 2.0 * etasq) - (2.0 * CK2) * tsi /
            (self.aodp * psisq) * (-3.0 * self.x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                                   0.75 * self.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
                                   np.cos(2.0 * self.omegao))))

        self.c5, self.c3, self.omgcof = 0.0, 0.0, 0.0

        if self.mode == SGDP4_NEAR_NORM:
            self.c5 = (2.0 * coef_1 * self.aodp * betao2 *
                       (1.0 + 2.75 * (etasq + eeta) + eeta * etasq))
            if self.eo > ECC_ALL:
                self.c3 = coef * tsi * A3OVK2 * \
                    self.xnodp * AE * self.sinIO / self.eo
            self.omgcof = self.bstar * self.c3 * np.cos(self.omegao)

        temp1 = 3.0 * CK2 * pinvsq * self.xnodp
        temp2 = temp1 * CK2 * pinvsq
        temp3 = 1.25 * CK4 * pinvsq**2 * self.xnodp

        self.xmdot = (self.xnodp + (0.5 * temp1 * betao * self.x3thm1 + 0.0625 *
                                    temp2 * betao * (13.0 - 78.0 * theta2 +
                                                     137.0 * theta4)))

        x1m5th = 1.0 - 5.0 * theta2

        self.omgdot = (-0.5 * temp1 * x1m5th + 0.0625 * temp2 *
                       (7.0 - 114.0 * theta2 + 395.0 * theta4) +
                       temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4))

        xhdot1 = -temp1 * self.cosIO
        self.xnodot = (xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) +
                                 2.0 * temp3 * (3.0 - 7.0 * theta2)) * self.cosIO)

        if self.eo > ECC_ALL:
            self.xmcof = (-(2. / 3) * AE) * coef * self.bstar / eeta
        else:
            self.xmcof = 0.0

        self.xnodcf = 3.5 * betao2 * xhdot1 * self.c1
        self.t2cof = 1.5 * self.c1

        # Check for possible divide-by-zero for X/(1+cos(xincl)) when
        # calculating xlcof */
        temp0 = 1.0 + self.cosIO
        if np.abs(temp0) < EPS_COS:
            temp0 = np.sign(temp0) * EPS_COS

        self.xlcof = 0.125 * A3OVK2 * self.sinIO * \
            (3.0 + 5.0 * self.cosIO) / temp0

        self.aycof = 0.25 * A3OVK2 * self.sinIO

        self.cosXMO = np.cos(self.xmo)
        self.sinXMO = np.sin(self.xmo)
        self.delmo = (1.0 + self.eta * self.cosXMO)**3

        if self.mode == SGDP4_NEAR_NORM:
            c1sq = self.c1**2
            self.d2 = 4.0 * self.aodp * tsi * c1sq
            temp0 = self.d2 * tsi * self.c1 / 3.0
            self.d3 = (17.0 * self.aodp + s4) * temp0
            self.d4 = 0.5 * temp0 * self.aodp * tsi * \
                (221.0 * self.aodp + 31.0 * s4) * self.c1
            self.t3cof = self.d2 + 2.0 * c1sq
            self.t4cof = 0.25 * \
                (3.0 * self.d3 + self.c1 * (12.0 * self.d2 + 10.0 * c1sq))
            self.t5cof = (0.2 * (3.0 * self.d4 + 12.0 * self.c1 * self.d3 + 6.0 * self.d2**2 +
                                 15.0 * c1sq * (2.0 * self.d2 + c1sq)))

        elif self.mode == SGDP4_DEEP_NORM:
            raise NotImplementedError('Deep space calculations not supported')

    def propagate(self, utc_time):
        kep = {}

        # get the time delta in minutes
        # ts = astronomy._days(utc_time - self.t_0) * XMNPDA
        # print utc_time.shape
        # print self.t_0
        utc_time = dt2np(utc_time)
        ts = (utc_time - self.t_0) / np.timedelta64(1, 'm')

        em = self.eo
        xinc = self.xincl

        xmp = self.xmo + self.xmdot * ts
        xnode = self.xnodeo + ts * (self.xnodot + ts * self.xnodcf)
        omega = self.omegao + self.omgdot * ts

        if self.mode == SGDP4_ZERO_ECC:
            raise NotImplementedError('Mode SGDP4_ZERO_ECC not implemented')
        elif self.mode == SGDP4_NEAR_SIMP:
            raise NotImplementedError('Mode "Near-space, simplified equations"'
                                      ' not implemented')
        elif self.mode == SGDP4_NEAR_NORM:
            delm = self.xmcof * \
                ((1.0 + self.eta * np.cos(xmp))**3 - self.delmo)
            temp0 = ts * self.omgcof + delm
            xmp += temp0
            omega -= temp0
            tempa = 1.0 - \
                (ts *
                 (self.c1 + ts * (self.d2 + ts * (self.d3 + ts * self.d4))))
            tempe = self.bstar * \
                (self.c4 * ts + self.c5 * (np.sin(xmp) - self.sinXMO))
            templ = ts * ts * \
                (self.t2cof + ts *
                 (self.t3cof + ts * (self.t4cof + ts * self.t5cof)))
            a = self.aodp * tempa**2
            e = em - tempe
            xl = xmp + omega + xnode + self.xnodp * templ

        else:
            raise NotImplementedError('Deep space calculations not supported')

        if np.any(a < 1):
            raise Exception('Satellite crashed at time %s', utc_time)
        elif np.any(e < ECC_LIMIT_LOW):
            raise ValueError('Satellite modified eccentricity too low: %s < %e'
                             % (str(e[e < ECC_LIMIT_LOW]), ECC_LIMIT_LOW))

        e = np.where(e < ECC_EPS, ECC_EPS, e)
        e = np.where(e > ECC_LIMIT_HIGH, ECC_LIMIT_HIGH, e)

        beta2 = 1.0 - e**2

        # Long period periodics
        sinOMG = np.sin(omega)
        cosOMG = np.cos(omega)

        temp0 = 1.0 / (a * beta2)
        axn = e * cosOMG
        ayn = e * sinOMG + temp0 * self.aycof
        xlt = xl + temp0 * self.xlcof * axn

        elsq = axn**2 + ayn**2

        if np.any(elsq >= 1):
            raise Exception('e**2 >= 1 at %s', utc_time)

        kep['ecc'] = np.sqrt(elsq)

        epw = np.fmod(xlt - xnode, 2 * np.pi)
        # needs a copy in case of an array
        capu = np.array(epw)
        maxnr = kep['ecc']
        for i in range(10):
            sinEPW = np.sin(epw)
            cosEPW = np.cos(epw)

            ecosE = axn * cosEPW + ayn * sinEPW
            esinE = axn * sinEPW - ayn * cosEPW
            f = capu - epw + esinE
            if np.all(np.abs(f) < NR_EPS):
                break

            df = 1.0 - ecosE

            # 1st order Newton-Raphson correction.
            nr = f / df

            # 2nd order Newton-Raphson correction.
            nr = np.where(np.logical_and(i == 0, np.abs(nr) > 1.25 * maxnr),
                          np.sign(nr) * maxnr,
                          f / (df + 0.5 * esinE * nr))
            epw += nr

        # Short period preliminary quantities
        temp0 = 1.0 - elsq
        betal = np.sqrt(temp0)
        pl = a * temp0
        r = a * (1.0 - ecosE)
        invR = 1.0 / r
        temp2 = a * invR
        temp3 = 1.0 / (1.0 + betal)
        cosu = temp2 * (cosEPW - axn + ayn * esinE * temp3)
        sinu = temp2 * (sinEPW - ayn - axn * esinE * temp3)
        u = np.arctan2(sinu, cosu)
        sin2u = 2.0 * sinu * cosu
        cos2u = 2.0 * cosu**2 - 1.0
        temp0 = 1.0 / pl
        temp1 = CK2 * temp0
        temp2 = temp1 * temp0

        # Update for short term periodics to position terms.

        rk = r * (1.0 - 1.5 * temp2 * betal * self.x3thm1) + \
            0.5 * temp1 * self.x1mth2 * cos2u
        uk = u - 0.25 * temp2 * self.x7thm1 * sin2u
        xnodek = xnode + 1.5 * temp2 * self.cosIO * sin2u
        xinck = xinc + 1.5 * temp2 * self.cosIO * self.sinIO * cos2u

        if np.any(rk < 1):
            raise Exception('Satellite crashed at time %s', utc_time)

        temp0 = np.sqrt(a)
        temp2 = XKE / (a * temp0)
        rdotk = ((XKE * temp0 * esinE * invR - temp2 * temp1 * self.x1mth2 * sin2u) *
                 (XKMPER / AE * XMNPDA / 86400.0))
        rfdotk = ((XKE * np.sqrt(pl) * invR + temp2 * temp1 *
                   (self.x1mth2 * cos2u + 1.5 * self.x3thm1)) *
                  (XKMPER / AE * XMNPDA / 86400.0))

        kep['radius'] = rk * XKMPER / AE
        kep['theta'] = uk
        kep['eqinc'] = xinck
        kep['ascn'] = xnodek
        kep['argp'] = omega
        kep['smjaxs'] = a * XKMPER / AE
        kep['rdotk'] = rdotk
        kep['rfdotk'] = rfdotk

        return kep


def kep2xyz(kep):
    sinT = np.sin(kep['theta'])
    cosT = np.cos(kep['theta'])
    sinI = np.sin(kep['eqinc'])
    cosI = np.cos(kep['eqinc'])
    sinS = np.sin(kep['ascn'])
    cosS = np.cos(kep['ascn'])

    xmx = -sinS * cosI
    xmy = cosS * cosI

    ux = xmx * sinT + cosS * cosT
    uy = xmy * sinT + sinS * cosT
    uz = sinI * sinT

    x = kep['radius'] * ux
    y = kep['radius'] * uy
    z = kep['radius'] * uz

    vx = xmx * cosT - cosS * sinT
    vy = xmy * cosT - sinS * sinT
    vz = sinI * cosT

    v_x = kep['rdotk'] * ux + kep['rfdotk'] * vx
    v_y = kep['rdotk'] * uy + kep['rfdotk'] * vy
    v_z = kep['rdotk'] * uz + kep['rfdotk'] * vz

    return np.array((x, y, z)), np.array((v_x, v_y, v_z))


if __name__ == "__main__":
    obs_lon, obs_lat = np.deg2rad((12.4143, 55.9065))
    obs_alt = 0.02
    o = Orbital(satellite="METOP-B")

    t_start = datetime.now()
    t_stop = t_start + timedelta(minutes=20)
    t = t_start
    while t < t_stop:
        t += timedelta(seconds=15)
        lon, lat, alt = o.get_lonlatalt(t)
        lon, lat = np.rad2deg((lon, lat))
        az, el = o.get_observer_look(t, obs_lon, obs_lat, obs_alt)
        ob = o.get_orbit_number(t, tbus_style=True)
        print(az, el, ob)
