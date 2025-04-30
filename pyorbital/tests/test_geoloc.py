#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2023 Pytroll Community

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

"""Test the geoloc module."""


import datetime as dt

import numpy as np
import pytest

from pyorbital.geoloc import ScanGeometry, geodetic_lat, qrotate, subpoint
from pyorbital.geoloc_avhrr import compute_avhrr_gcps_lonlatalt, estimate_time_and_attitude_deviations
from pyorbital.geoloc_instrument_definitions import (
    amsua,
    ascat,
    atms,
    avhrr,
    avhrr_from_times,
    avhrr_gac_from_times,
    hirs4,
    mhs,
    slstr_nadir,
    viirs,
)


class TestQuaternion:
    """Test the quaternion rotation."""

    def test_qrotate(self):
        """Test quaternion rotation."""
        vector = np.array([[1, 0, 0]]).T
        axis = np.array([[0, 1, 0]]).T
        angle = np.deg2rad(90)

        result = qrotate(vector, axis, angle)[:, 0]
        expected = np.array([0, 0, 1])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

        axis = np.array([0, 1, 0])
        result = qrotate(vector, axis, angle)
        expected = np.array([[0, 0, 1]]).T
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

        vector = np.array([[1, 0, 0],
                           [0, 0, 1]]).T
        axis = np.array([0, 1, 0])
        angle = np.deg2rad(90)
        result = qrotate(vector, axis, angle)
        expected = np.array([[0, 0, 1],
                             [-1, 0, 0]]).T

        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

        axis = np.array([[0, 1, 0]]).T
        result = qrotate(vector, axis, angle)
        expected = np.array([[0, 0, 1],
                             [-1, 0, 0]]).T

        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)


class TestGeoloc:
    """Test for the core computing part."""

    def test_scan_geometry(self):
        """Test the ScanGeometry object."""
        scans_nb = 1

        xy = np.vstack((np.deg2rad(np.array([10, 0, -10])),
                        np.array([0, 0, 0])))
        xy = np.tile(xy[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

        times = np.tile([-0.1, 0, 0.1], [np.int32(scans_nb), 1])

        instrument = ScanGeometry(xy, times)

        np.testing.assert_allclose(np.rad2deg(instrument.fovs[0]), np.array([[10, 0, -10]]))

        # Test vectors

        pos = np.rollaxis(np.tile(np.array([0, 0, 7000]), [3, 1, 1]), 2)
        vel = np.rollaxis(np.tile(np.array([1, 0, 0]), [3, 1, 1]), 2)
        pos = np.stack([np.array([0, 0, 7000])] * 3, 1)[:, np.newaxis, :]
        vel = np.stack([np.array([1, 0, 0])] * 3, 1)[:, np.newaxis, :]

        vec = instrument.vectors(pos, vel)

        result = vec[:, 0, 1]
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

        # Check if we can pass an array for yaw
        vec = instrument.vectors(pos, vel, yaw=[0])

        result = vec[:, 0, 1]
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)


        # minus sin because we use trigonometrical direction of angles
        result = vec[:, 0, 0]
        expected = np.array([0, -np.sin(np.deg2rad(10)), -np.cos(np.deg2rad(10))])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)

        result = vec[:, 0, 2]
        expected = np.array([0, -np.sin(np.deg2rad(-10)), -np.cos(np.deg2rad(-10))])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)

        # Test times

        start_of_scan = np.datetime64(dt.datetime(2014, 1, 8, 11, 30))
        times = instrument.times(start_of_scan)

        assert times[0, 1] == start_of_scan
        assert times[0, 0] == start_of_scan - np.timedelta64(100, "ms")
        assert times[0, 2] == start_of_scan + np.timedelta64(100, "ms")

    def test_geodetic_lat(self):
        """Test the determination of the geodetic latitude."""
        point = np.array([[7000, 0, 7000]]).T
        np.testing.assert_allclose(geodetic_lat(point),
                                   np.array([0.78755832699854733]), rtol=1e-8, atol=1e-8)

        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        result = geodetic_lat(points)
        expected = np.array([0.78755832699854733, 0.78755832699854733])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

    def test_subpoint(self):
        """Test nadir determination."""
        a = 6378.137  # km
        b = 6356.75231414  # km, GRS80
        point = np.array([0, 0, 7000])
        nadir = subpoint(point, a, b)
        np.testing.assert_allclose(nadir, np.array([0, 0, b]), rtol=1e-7, atol=1e-7)

        point = np.array([7000, 0, 7000])
        nadir = subpoint(point, a, b)
        np.testing.assert_allclose(nadir,
                                   np.array([4507.85431429,
                                             0,
                                             4497.06396339]), rtol=1e-8, atol=1e-8)
        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        nadir = subpoint(points, a, b)
        np.testing.assert_allclose(nadir[:, 0],
                                   np.array([4507.85431429,
                                             0,
                                             4497.06396339]), rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(nadir[:, 1],
                                   np.array([4507.85431429,
                                             0,
                                             4497.06396339]), rtol=1e-8, atol=1e-8)




def test_arbitrary_point_geoloc():
    """Test geolocating an arbitrary point in the swath."""
    from pyorbital.geoloc_avhrr import compute_avhrr_gcps_lonlatalt

    # Couple of example Two Line Elements
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"

    # Choosing a specific time, this should be relatively close to the issue date of the TLE
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)
    rpy = (0, 0, 0)

    max_scan_angle = 55.37

    gcps = np.array([[2, 500], [1500, 700], [20, 1000]])

    lons, lats, alts = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, t, (tle1, tle2))

    assert lons[0] == pytest.approx(-34.69996894)
    assert lats[0] == pytest.approx(56.69799502)

    assert lons[2] == pytest.approx(-27.573052737698944)
    assert lats[2] == pytest.approx(55.626740897592654)


def test_minimize_geoloc_error():
    """Test minimizing the distance to a set of gcps."""
    # Couple of example Two Line Elements
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    tle = (tle1, tle2)

    # Choosing a specific time, this should be relatively close to the issue date of the TLE
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    ref_time_displacement = 0.51
    ref_time = t + dt.timedelta(seconds=ref_time_displacement)
    ref_yaw = 0.1
    rpy = (0, 0, ref_yaw)
    max_scan_angle = 55.37
    # gcps are line/col
    gcps = np.array([[2, 500], [1500, 700], [20, 1000], [500, 1100], [100, 2000]])
    ref_lons, ref_lats, _ = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, ref_time, tle)
    time_diff, roll, pitch, yaw = estimate_time_and_attitude_deviations(gcps, ref_lons, ref_lats, t,
                                                                        tle, max_scan_angle)
    assert time_diff == pytest.approx(ref_time_displacement, abs=1e-2)
    assert yaw == pytest.approx(ref_yaw, abs=1e-2)


class TestGeolocDefs:
    """Test the instrument definitions."""

    def test_avhrr(self):
        """Test the definition of the avhrr instrument."""
        avh = avhrr(1, np.array([0, 1023.5, 2047]))
        result = np.rad2deg(avh.fovs[0])
        expected = np.array([[55.37, 0, -55.37]])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)

        avh = avhrr(1, np.array([0, 1023.5, 2047]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[10, 0, -10]]))

        # This is perhaps a bit odd, to require avhrr to accept floats for
        # the number of scans? FIXME!
        avh = avhrr(1.1, np.array([0, 1023.5, 2047]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[10, 0, -10]]))

    def test_avhrr_from_times(self):
        """Test generating the avhrr from times."""
        avh = avhrr_from_times([dt.datetime(2000,1,1,0,0,0)], [0, 1023.5, 2047])
        result = np.rad2deg(avh.fovs[0])
        expected = np.array([[55.37, 0, -55.37]])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)
        result = avh.times(dt.date(2000,1,1))
        expected = ((np.array([[0, 1023.5, 2047]]) * 25000).astype("timedelta64[ns]")
                    + np.datetime64("2000-01-01T00:00:00"))
        np.testing.assert_equal(result, expected)

        avh = avhrr_from_times([dt.datetime(2000,1,1,0,0,0)], np.array([0, 1023.5, 2047]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[10, 0, -10]]))

        avh = avhrr_from_times([dt.datetime(2000,1,1,0,0,0), dt.datetime(2000,1,1,0,1,0)],
                    [0, 2047])
        times = avh.times(dt.datetime(2001,1,1))
        expected = (np.array([[0,51175000],[60000000000, 60051175000]]).astype("timedelta64[ns]")
                    + np.datetime64("2001-01-01"))
        np.testing.assert_equal(times, expected)


    def test_avhrr_gac_from_times(self):
        """Test getting avhrr gac from times."""
        avh = avhrr_gac_from_times([dt.datetime(2000,1,1,0,0,0)], [0, 204, 408])
        result = np.rad2deg(avh.fovs[0])
        expected = np.array([[55.180655, 0, -55.180655]])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)
        result = avh.times(dt.date(2000,1,1))
        expected = ((np.array([[0, 204, 408]]) * 125000).astype("timedelta64[ns]")
                    + np.datetime64("2000-01-01T00:00:00"))
        np.testing.assert_equal(result, expected)

        avh = avhrr_gac_from_times([dt.datetime(2000,1,1,0,0,0)], np.array([0, 204, 408]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[9.965804, 0, -9.965804]]))

        avh = avhrr_gac_from_times([dt.datetime(2000,1,1,0,0,0), dt.datetime(2000,1,1,0,1,0)],
                    [0, 408])
        times = avh.times(dt.datetime(2001,1,1))
        expected = (np.array([[0,51000000],[60000000000, 60051000000]]).astype("timedelta64[ns]")
                    + np.datetime64("2001-01-01"))
        np.testing.assert_equal(times, expected)


    def test_viirs(self):
        """Test the definition of the viirs instrument."""
        geom = viirs(1, np.array([0, 3200, 6399]))
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32, 1]),
            np.tile(np.array([[0., -0., 0]]), [32, 1])], dtype=np.float64)

        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

        geom = viirs(2, np.array([0, 3200, 6399]))
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32*2, 1]),
            np.tile(np.array([[0., -0., 0]]), [32*2, 1])], dtype=np.float64)

        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_viirs_defaults(self):
        """Test the definition of the viirs instrument with default slicing."""
        geom = viirs(1, chn_pixels=3)
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32, 1]),
            np.tile(np.array([[0., -0., 0]]), [32, 1])], dtype=np.float64)

        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_amsua(self):
        """Test the definition of the amsua instrument."""
        geom = amsua(1)
        expected_fovs = np.array([
            [[0.84,  0.78,  0.73,  0.67,  0.61,  0.55,  0.49,  0.44,  0.38,
              0.32,  0.26,  0.2,  0.15,  0.09,  0.03, -0.03, -0.09, -0.15,
              -0.2, -0.26, -0.32, -0.38, -0.44, -0.49, -0.55, -0.61, -0.67,
              -0.73, -0.78, -0.84]],
            np.zeros((1, 30))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

    def test_mhs(self):
        """Test the definition of the mhs instrument."""
        geom = mhs(1)
        expected_fovs = np.array([
            [[0.86,  0.84,  0.82,  0.8,  0.79,  0.77,  0.75,  0.73,  0.71,
              0.69,  0.67,  0.65,  0.63,  0.61,  0.59,  0.57,  0.55,  0.53,
              0.51,  0.49,  0.48,  0.46,  0.44,  0.42,  0.4,  0.38,  0.36,
              0.34,  0.32,  0.3,  0.28,  0.26,  0.24,  0.22,  0.2,  0.18,
              0.16,  0.15,  0.13,  0.11,  0.09,  0.07,  0.05,  0.03,  0.01,
              -0.01, -0.03, -0.05, -0.07, -0.09, -0.11, -0.13, -0.15, -0.16,
              -0.18, -0.2, -0.22, -0.24, -0.26, -0.28, -0.3, -0.32, -0.34,
              -0.36, -0.38, -0.4, -0.42, -0.44, -0.46, -0.48, -0.49, -0.51,
              -0.53, -0.55, -0.57, -0.59, -0.61, -0.63, -0.65, -0.67, -0.69,
              -0.71, -0.73, -0.75, -0.77, -0.79, -0.8, -0.82, -0.84, -0.86]],
            np.zeros((1, 90))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_hirs4(self):
        """Test the definition of the hirs4 instrument."""
        geom = hirs4(1)
        expected_fovs = np.array([
            [[0.86,  0.83,  0.8,  0.77,  0.74,  0.71,  0.68,  0.64,  0.61,
              0.58,  0.55,  0.52,  0.49,  0.46,  0.42,  0.39,  0.36,  0.33,
              0.3,  0.27,  0.24,  0.2,  0.17,  0.14,  0.11,  0.08,  0.05,
              0.02, -0.02, -0.05, -0.08, -0.11, -0.14, -0.17, -0.2, -0.24,
              -0.27, -0.3, -0.33, -0.36, -0.39, -0.42, -0.46, -0.49, -0.52,
              -0.55, -0.58, -0.61, -0.64, -0.68, -0.71, -0.74, -0.77, -0.8,
              -0.83, -0.86]],
            np.zeros((1, 56))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_atms(self):
        """Test the definition of the atms instrument."""
        geom = atms(1)
        expected_fovs = np.array([
            [[0.92,  0.9,  0.88,  0.86,  0.84,  0.82,  0.8,  0.78,  0.76,
              0.75,  0.73,  0.71,  0.69,  0.67,  0.65,  0.63,  0.61,  0.59,
              0.57,  0.55,  0.53,  0.51,  0.49,  0.47,  0.46,  0.44,  0.42,
              0.4,  0.38,  0.36,  0.34,  0.32,  0.3,  0.28,  0.26,  0.24,
              0.22,  0.2,  0.18,  0.16,  0.15,  0.13,  0.11,  0.09,  0.07,
              0.05,  0.03,  0.01, -0.01, -0.03, -0.05, -0.07, -0.09, -0.11,
              -0.13, -0.15, -0.16, -0.18, -0.2, -0.22, -0.24, -0.26, -0.28,
              -0.3, -0.32, -0.34, -0.36, -0.38, -0.4, -0.42, -0.44, -0.46,
              -0.47, -0.49, -0.51, -0.53, -0.55, -0.57, -0.59, -0.61, -0.63,
              -0.65, -0.67, -0.69, -0.71, -0.73, -0.75, -0.76, -0.78, -0.8,
              -0.82, -0.84, -0.86, -0.88, -0.9, -0.92]],
            np.zeros((1, 96))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_ascat(self):
        """Test the definition of the ASCAT instrument onboard Metop."""
        geom = ascat(1)
        expected_fovs = np.array([
            [[0.9250245,  0.90058989,  0.87615528,  0.85172067,
              0.82728607,  0.80285146,  0.77841685,  0.75398224,
              0.72954763,  0.70511302,  0.68067841,  0.6562438,
              0.63180919,  0.60737458,  0.58293997,  0.55850536,
              0.53407075,  0.50963614,  0.48520153,  0.46076692,
              0.43633231, -0.43633231, -0.46076692, -0.48520153,
              -0.50963614, -0.53407075, -0.55850536, -0.58293997,
              -0.60737458, -0.63180919, -0.6562438, -0.68067841,
              -0.70511302, -0.72954763, -0.75398224, -0.77841685,
              -0.80285146, -0.82728607, -0.85172067, -0.87615528,
              -0.90058989, -0.9250245]], np.zeros((1, 42))], dtype=np.float64)

        np.testing.assert_allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)
        geom = ascat(1, np.array([0, 41]))
        expected_fovs = np.array([[[0.9250245,  -0.9250245]],
                                  [[0.,  0.]]], dtype=np.float64)
        np.testing.assert_allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

        geom = ascat(1, np.array([0, -1]))
        np.testing.assert_allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

    def test_slstr_nadir(self):
        """Test the definition of the slstr instrument nadir view flying on Sentinel-3."""
        geom = slstr_nadir(1, [0, 1])

        expected_fovs = np.array([
            np.tile(np.array([[0.8115781, -0.38571776]]), [1, 1]),
            np.tile(np.array([[0., 0.]]), [1, 1])], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

        geom = slstr_nadir(1, None)

        np.testing.assert_equal(geom.fovs.size, 6000)
