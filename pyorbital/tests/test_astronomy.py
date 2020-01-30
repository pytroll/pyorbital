#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

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

import unittest

from datetime import datetime
import pyorbital.astronomy as astr


class TestAstronomy(unittest.TestCase):

    def setUp(self):
        pass

    def test_jdays(self):
        """Test julian day functions."""
        t = datetime(2000, 1, 1, 12, 0)
        self.assertEqual(astr.jdays(t), 2451545.0)
        self.assertEqual(astr.jdays2000(t), 0)
        t = datetime(2009, 10, 8, 14, 30)
        self.assertEqual(astr.jdays(t), 2455113.1041666665)
        self.assertEqual(astr.jdays2000(t), 3568.1041666666665)

    def test_sunangles(self):
        """Test the sun-angle calculations."""
        lat, lon = 58.6167, 16.1833  # Norrkoping
        time_slot = datetime(2011, 9, 23, 12, 0)

        sun_theta = astr.sun_zenith_angle(time_slot, lon, lat)
        self.assertAlmostEqual(sun_theta, 60.371433482557833, places=8)
        sun_theta = astr.sun_zenith_angle(time_slot, 0., 0.)
        self.assertAlmostEqual(sun_theta, 1.8751916863323426, places=8)


def suite():
    """The suite for test_astronomy."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestAstronomy))

    return mysuite
