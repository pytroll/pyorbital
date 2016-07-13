#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Martin Raspaud

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
from pyorbital.geoloc_instrument_definitions import avhrr
from pyorbital.geoloc import (ScanGeometry, qrotate,
                              subpoint, geodetic_lat)
import numpy as np
from datetime import datetime, timedelta

class TestQuaternion(unittest.TestCase):
    """Test the quaternion rotation.
    """

    def test_qrotate(self):
        """Test quaternion rotation
        """
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
        instrument = ScanGeometry(np.deg2rad(np.array([[10, 0],
                                                       [0, 0],
                                                       [-10, 0]])),
                                  np.array([-0.1, 0, 0.1]))
        self.assertTrue(np.allclose(np.rad2deg(instrument.fovs[:, 0]),
                                    np.array([10, 0, -10])))

        # Test vectors
        
        pos = np.array([[0, 0, 7000]]).T
        vel = np.array([[1, 0, 0]]).T
        vec = instrument.vectors(pos, vel)
        self.assertTrue(np.allclose(np.array([[0, 0, -1]]),
                                    vec[:, 1]))

        # minus sin because we use trigonometrical direction of angles
        
        self.assertTrue(np.allclose(np.array([[0,
                                               -np.sin(np.deg2rad(10)),
                                               -np.cos(np.deg2rad(10))]]),
                                    vec[:, 0]))
        self.assertTrue(np.allclose(np.array([[0,
                                               -np.sin(np.deg2rad(-10)),
                                               -np.cos(np.deg2rad(-10))]]),
                                    vec[:, 2]))

        # Test times

        start_of_scan = datetime(2014, 1, 8, 11, 30)
        times = instrument.times(start_of_scan)
        self.assertEquals(times[1], start_of_scan)
        self.assertEquals(times[0], start_of_scan - timedelta(seconds=0.1))
        self.assertEquals(times[2], start_of_scan + timedelta(seconds=0.1))

    def test_geodetic_lat(self):
        """Test the determination of the geodetic latitude.
        """

        a = 6378.137 # km
        b = 6356.75231414 # km, GRS80

        point = np.array([7000, 0, 7000])
        self.assertEqual(geodetic_lat(point), 0.78755832699854733)
        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        self.assertTrue(np.allclose(geodetic_lat(points),
                                    np.array([0.78755832699854733,
                                              0.78755832699854733])))
        
        
    def test_subpoint(self):
        """Test nadir determination.
        """
        a = 6378.137 # km
        b = 6356.75231414 # km, GRS80
        point = np.array([0, 0, 7000])
        nadir = subpoint(point, a, b)
        self.assertTrue(np.allclose(nadir, np.array([[0, 0, b]]).T))

        point = np.array([7000, 0, 7000])
        nadir = subpoint(point, a, b)
        self.assertTrue(np.allclose(nadir,
                                    np.array([[4507.85431429,
                                               0,
                                               4497.06396339]]).T))
        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        nadir = subpoint(points, a, b)
        self.assertTrue(np.allclose(nadir,
                                    np.array([[4507.85431429,
                                               0,
                                               4497.06396339]]).T))


        

class TestGeolocDefs(unittest.TestCase):
    """Test the instrument definitions.
    """
    def test_avhrr(self):
        """Test the definition of the avhrr instrument
        """
        avh = avhrr(1, np.array([0, 1023.5 ,2047]))
        self.assertTrue(np.allclose(np.rad2deg(avh.fovs[:, 0]),
                                    np.array([55.37, 0, -55.37])))

        avh = avhrr(1, np.array([0, 1023.5 ,2047]), 10)
        self.assertTrue(np.allclose(np.rad2deg(avh.fovs[:, 0]),
                                    np.array([10, 0, -10])))

def suite():
    """The suite for test_geoloc
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestQuaternion))
    mysuite.addTest(loader.loadTestsFromTestCase(TestGeoloc))
    mysuite.addTest(loader.loadTestsFromTestCase(TestGeolocDefs))
    
    return mysuite
