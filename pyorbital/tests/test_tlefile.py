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

"""Testing TLE file reading
"""


from pyorbital.tlefile import Tle
import datetime
import unittest

line0 = "ISS (ZARYA)"
line1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"


class TLETest(unittest.TestCase):
    """Test TLE reading.

    We're using the wikipedia example::

     ISS (ZARYA)
     1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
     2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

    """

    def check_example(self, tle):
        """Check the *tle* instance against predetermined values.
        """
        # line 1
        self.assertEqual(tle.satnumber, "25544")
        self.assertEqual(tle.classification, "U")
        self.assertEqual(tle.id_launch_year, "98")
        self.assertEqual(tle.id_launch_number, "067")
        self.assertEqual(tle.id_launch_piece.strip(), "A")
        self.assertEqual(tle.epoch_year, "08")
        self.assertEqual(tle.epoch_day, 264.51782528)
        epoch = (datetime.datetime(2008, 1, 1)
                 + datetime.timedelta(days=264.51782528 - 1))
        self.assertEqual(tle.epoch, epoch)
        self.assertEqual(tle.mean_motion_derivative, -.00002182)
        self.assertEqual(tle.mean_motion_sec_derivative, 0.0)
        self.assertEqual(tle.bstar, -.11606e-4)
        self.assertEqual(tle.ephemeris_type, 0)
        self.assertEqual(tle.element_number, 292)

        # line 2
        self.assertEqual(tle.inclination, 51.6416)
        self.assertEqual(tle.right_ascension, 247.4627)
        self.assertEqual(tle.excentricity, .0006703)
        self.assertEqual(tle.arg_perigee, 130.5360)
        self.assertEqual(tle.mean_anomaly, 325.0288)
        self.assertEqual(tle.mean_motion, 15.72125391)
        self.assertEqual(tle.orbit, 56353)

    def test_from_line(self):
        tle = Tle("ISS (ZARYA)", line1=line1, line2=line2)
        self.check_example(tle)

    def test_from_file(self):
        from tempfile import mkstemp
        from os import write, close, remove
        filehandle, filename = mkstemp()
        try:
            write(filehandle, "\n".join([line0, line1, line2]).encode('utf-8'))
            close(filehandle)
            tle = Tle("ISS (ZARYA)", filename)
            self.check_example(tle)
        finally:
            remove(filename)


def suite():
    """The suite for test_tlefile
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TLETest))

    return mysuite
