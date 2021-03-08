#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2017, 2018, 2021 Martin Raspaud

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
"""Test the geoloc module.
"""

import unittest
from datetime import datetime

import numpy as np

from pyorbital.geoloc import ScanGeometry, geodetic_lat, qrotate, subpoint
from pyorbital.geoloc_instrument_definitions import avhrr, viirs, amsua, mhs, hirs4, atms, ascat


class TestQuaternion(unittest.TestCase):
    """Test the quaternion rotation."""

    def test_qrotate(self):
        """Test quaternion rotation."""
        vector = np.array([[1, 0, 0]]).T
        axis = np.array([[0, 1, 0]]).T
        angle = np.deg2rad(90)
        self.assertTrue(np.allclose(qrotate(vector, axis, angle),
                                    np.array([[0, 0, 1]]).T))

        axis = np.array([0, 1, 0])
        self.assertTrue(np.allclose(qrotate(vector, axis, angle),
                                    np.array([[0, 0, 1]]).T))

        vector = np.array([[1, 0, 0],
                           [0, 0, 1]]).T
        axis = np.array([0, 1, 0])
        angle = np.deg2rad(90)
        self.assertTrue(np.allclose(qrotate(vector, axis, angle),
                                    np.array([[0, 0, 1],
                                              [-1, 0, 0]]).T))

        axis = np.array([[0, 1, 0]]).T
        self.assertTrue(np.allclose(qrotate(vector, axis, angle),
                                    np.array([[0, 0, 1],
                                              [-1, 0, 0]]).T))


class TestGeoloc(unittest.TestCase):

    """Test for the core computing part.
    """

    def test_scan_geometry(self):
        """Test the ScanGeometry object.
        """
        scans_nb = 1

        xy = np.vstack((np.deg2rad(np.array([10, 0, -10])),
                        np.array([0, 0, 0])))
        xy = np.tile(xy[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

        times = np.tile([-0.1, 0, 0.1], [np.int32(scans_nb), 1])

        instrument = ScanGeometry(xy, times)

        self.assertTrue(np.allclose(np.rad2deg(instrument.fovs[0]),
                                    np.array([[10, 0, -10]])))

        # Test vectors

        pos = np.rollaxis(np.tile(np.array([0, 0, 7000]), [3, 1, 1]), 2)
        vel = np.rollaxis(np.tile(np.array([1, 0, 0]), [3, 1, 1]), 2)
        pos = np.stack([np.array([0, 0, 7000])] * 3, 1)[:, np.newaxis, :]
        vel = np.stack([np.array([1, 0, 0])] * 3, 1)[:, np.newaxis, :]

        vec = instrument.vectors(pos, vel)

        self.assertTrue(np.allclose(np.array([[0, 0, -1]]),
                                    vec[:, 0, 1]))

        # minus sin because we use trigonometrical direction of angles

        self.assertTrue(np.allclose(np.array([[0,
                                               -np.sin(np.deg2rad(10)),
                                               -np.cos(np.deg2rad(10))]]),
                                    vec[:, 0, 0]))
        self.assertTrue(np.allclose(np.array([[0,
                                               -np.sin(np.deg2rad(-10)),
                                               -np.cos(np.deg2rad(-10))]]),
                                    vec[:, 0, 2]))

        # Test times

        start_of_scan = np.datetime64(datetime(2014, 1, 8, 11, 30))
        times = instrument.times(start_of_scan)

        self.assertEqual(times[0, 1], start_of_scan)
        self.assertEqual(times[0, 0], start_of_scan -
                         np.timedelta64(100, 'ms'))
        self.assertEqual(times[0, 2], start_of_scan +
                         np.timedelta64(100, 'ms'))

    def test_geodetic_lat(self):
        """Test the determination of the geodetic latitude."""
        point = np.array([7000, 0, 7000])
        self.assertEqual(geodetic_lat(point), 0.78755832699854733)
        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        self.assertTrue(np.allclose(geodetic_lat(points), np.array([0.78755832699854733, 0.78755832699854733])))

    def test_subpoint(self):
        """Test nadir determination."""
        a = 6378.137  # km
        b = 6356.75231414  # km, GRS80
        point = np.array([0, 0, 7000])
        nadir = subpoint(point, a, b)
        self.assertTrue(np.allclose(nadir, np.array([[0, 0, b]])))

        point = np.array([7000, 0, 7000])
        nadir = subpoint(point, a, b)
        self.assertTrue(np.allclose(nadir,
                                    np.array([[4507.85431429,
                                               0,
                                               4497.06396339]])))
        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        nadir = subpoint(points, a, b)
        self.assertTrue(np.allclose(nadir[:, 0],
                                    np.array([[4507.85431429,
                                               0,
                                               4497.06396339]])))
        self.assertTrue(np.allclose(nadir[:, 1],
                                    np.array([[4507.85431429,
                                               0,
                                               4497.06396339]])))


class TestGeolocDefs(unittest.TestCase):

    """Test the instrument definitions.
    """

    def test_avhrr(self):
        """Test the definition of the avhrr instrument
        """
        avh = avhrr(1, np.array([0, 1023.5, 2047]))
        self.assertTrue(np.allclose(np.rad2deg(avh.fovs[0]),
                                    np.array([55.37, 0, -55.37])))

        avh = avhrr(1, np.array([0, 1023.5, 2047]), 10)
        self.assertTrue(np.allclose(np.rad2deg(avh.fovs[0]),
                                    np.array([10, 0, -10])))

        # This is perhaps a bit odd, to require avhrr to accept floats for
        # the number of scans? FIXME!
        avh = avhrr(1.1, np.array([0, 1023.5, 2047]), 10)
        self.assertTrue(np.allclose(np.rad2deg(avh.fovs[0]),
                                    np.array([10, 0, -10])))

    def test_viirs(self):
        """Test the definition of the viirs instrument
        """
        geom = viirs(1, np.array([0, 3200, 6399]))
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32, 1]),
            np.tile(np.array([[0., -0., 0]]), [32, 1])], dtype=np.float64)

        self.assertTrue(np.allclose(geom.fovs,
                                    expected_fovs, rtol=1e-2, atol=1e-2))

        geom = viirs(2, np.array([0, 3200, 6399]))
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32*2, 1]),
            np.tile(np.array([[0., -0., 0]]), [32*2, 1])], dtype=np.float64)

        self.assertTrue(np.allclose(geom.fovs,
                                    expected_fovs, rtol=1e-2, atol=1e-2))

    def test_amsua(self):
        """Test the definition of the amsua instrument
        """
        geom = amsua(1)
        expected_fovs = np.array([
            [[0.84,  0.78,  0.73,  0.67,  0.61,  0.55,  0.49,  0.44,  0.38,
              0.32,  0.26,  0.2,  0.15,  0.09,  0.03, -0.03, -0.09, -0.15,
              -0.2, -0.26, -0.32, -0.38, -0.44, -0.49, -0.55, -0.61, -0.67,
              -0.73, -0.78, -0.84]],
            np.zeros((1, 30))], dtype=np.float64)
        self.assertTrue(np.allclose(geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2))

    def test_mhs(self):
        """Test the definition of the mhs instrument
        """
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
        self.assertTrue(np.allclose(geom.fovs,
                                    expected_fovs, rtol=1e-2, atol=1e-2))

    def test_hirs4(self):
        """Test the definition of the hirs4 instrument
        """
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
        self.assertTrue(np.allclose(geom.fovs,
                                    expected_fovs, rtol=1e-2, atol=1e-2))

    def test_atms(self):
        """Test the definition of the atms instrument
        """
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
        self.assertTrue(np.allclose(geom.fovs,
                                    expected_fovs, rtol=1e-2, atol=1e-2))

    def test_ascat(self):
        """Test the definition of the ASCAT instrument onboard Metop"""

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

        self.assertTrue(np.allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2))
        geom = ascat(1, np.array([0, 41]))
        expected_fovs = np.array([[[0.9250245,  -0.9250245]],
                                  [[0.,  0.]]], dtype=np.float64)
        self.assertTrue(np.allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2))

        geom = ascat(1, np.array([0, -1]))
        self.assertTrue(np.allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2))


def suite():
    """The suite for test_geoloc
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestQuaternion))
    mysuite.addTest(loader.loadTestsFromTestCase(TestGeoloc))
    mysuite.addTest(loader.loadTestsFromTestCase(TestGeolocDefs))

    return mysuite
