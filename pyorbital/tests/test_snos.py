#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c21529.ad.smhi.se>

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

"""Test the functions to find SNOs from two different satellite platforms."""

from pyorbital.orbital import Orbital
from pyorbital.sno_utils import get_arc
from pyorbital.sno_utils import get_sno_point
import sys
from datetime import datetime, timedelta

if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest


class GetSnosTest(unittest.TestCase):
    """Test getting snos."""

    def setUp(self):
        """Set up the test."""
        self.platform_id = 'Metop-B'
        self.orb_ref1 = Orbital('EOS-Aqua',
                                line1='1 27424U 02022A   14365.22903259  .00000847  00000-0  19796-3 0  2320',
                                line2='2 27424  98.2234 303.5330 0001197  87.2691  48.8815 14.57121895673369')
        self.orb_cmp1 = Orbital('Metop-B',
                                line1='1 38771U 12049A   14364.85213337  .00000110  00000-0  70147-4 0  6893',
                                line2='2 38771  98.7220  62.0425 0000591  33.5587  96.4312 14.21479121118503')
        self.tobj1 = datetime(2015, 1, 1, 0, 20)
        self.delta_t = timedelta(seconds=1200)
        self.timestep = timedelta(seconds=600)
        self.minthr = 10
        self.station = {'lon': 16.1465, 'lat': 58.578, 'alt': 0.03}

        self.orb_ref2 = Orbital('EOS-Aqua',
                                line1='1 27424U 02022A   15003.03711450  .00000723  00000-0  17040-3 0  2344',
                                line2='2 27424  98.2239 306.3049 0001165  87.8737   9.7700 14.57125023673775')
        self.orb_cmp2 = Orbital('Metop-B',
                                line1='1 38771U 12049A   15002.73354670  .00000116  00000-0  73308-4 0  6917',
                                line2='2 38771  98.7211  64.8896 0000438  20.9656  85.8395 14.21479960118911')
        self.tobj2 = datetime(2015, 1, 3, 17, 20)

    def test_get_sno_point(self):
        """Test getting a SNO point."""
        arc_ref_pltfrm = get_arc(self.tobj1, self.timestep, self.orb_ref1)
        arc_cmp_pltfrm = get_arc(self.tobj1, self.timestep, self.orb_cmp1)
        sno = get_sno_point(self.platform_id, self.orb_cmp1, self.orb_ref1,
                            self.tobj1, self.delta_t, arc_ref_pltfrm,
                            arc_cmp_pltfrm, self.station, self.minthr)

        expected = {'maxt': datetime(2015, 1, 1, 0, 12, 0, 506020),
                    'maxt_ref_pltfrm': datetime(2015, 1, 1, 0, 10, 36, 45923),
                    'ref_pltfrmsec': 36.045923,
                    'sec': 0.50602,
                    'is_within_antenna_horizon': False,
                    'tdmin': -1.4166666666666667,
                    'orbit_nr': 11866,
                    'point': (-9.057282984227129, -74.5159273899483)}

        self.assertAlmostEqual(expected['sec'], sno['sec'])
        self.assertAlmostEqual(expected['ref_pltfrmsec'], sno['ref_pltfrmsec'])
        self.assertAlmostEqual(expected['tdmin'], sno['tdmin'])
        self.assertTrue(expected['is_within_antenna_horizon'] == sno['is_within_antenna_horizon'])
        self.assertEqual(expected['orbit_nr'], sno['orbit_nr'])
        self.assertEqual(expected['maxt'], sno['maxt'])
        self.assertEqual(expected['maxt_ref_pltfrm'], sno['maxt_ref_pltfrm'])
        self.assertAlmostEqual(expected['point'][0], sno['point'][0])
        self.assertAlmostEqual(expected['point'][1], sno['point'][1])

        arc_ref_pltfrm = get_arc(self.tobj2, self.timestep, self.orb_ref2)
        arc_cmp_pltfrm = get_arc(self.tobj2, self.timestep, self.orb_cmp2)
        sno = get_sno_point(self.platform_id, self.orb_cmp2, self.orb_ref2,
                            self.tobj2, self.delta_t, arc_ref_pltfrm,
                            arc_cmp_pltfrm, self.station, self.minthr)

        expected = {'maxt': datetime(2015, 1, 3, 17, 13, 54, 493133),
                    'maxt_ref_pltfrm': datetime(2015, 1, 3, 17, 16, 46, 332581),
                    'ref_pltfrmsec': 46.332581,
                    'sec': 54.493133,
                    'is_within_antenna_horizon': False,
                    'tdmin': 2.85,
                    'orbit_nr': 11905,
                    'point': (-82.65813724524172, 76.0882567980114)}

        self.assertAlmostEqual(expected['sec'], sno['sec'])
        self.assertAlmostEqual(expected['ref_pltfrmsec'], sno['ref_pltfrmsec'])
        self.assertAlmostEqual(expected['tdmin'], sno['tdmin'])
        self.assertTrue(expected['is_within_antenna_horizon'] == sno['is_within_antenna_horizon'])
        self.assertEqual(expected['orbit_nr'], sno['orbit_nr'])
        self.assertEqual(expected['maxt'], sno['maxt'])
        self.assertEqual(expected['maxt_ref_pltfrm'], sno['maxt_ref_pltfrm'])
        self.assertAlmostEqual(expected['point'][0], sno['point'][0])
        self.assertAlmostEqual(expected['point'][1], sno['point'][1])


def suite():
    """Run the test suite for testing finding SNOs."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(GetSnosTest))

    return mysuite
