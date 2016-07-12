import unittest
from datetime import datetime

import numpy as np

from pyorbital import orbital

eps_deg = 10e-3


class TestOrbital(unittest.TestCase):

    def test_get_orbit_number(self):
        """Testing getting the orbitnumber from the tle"""
        sat = orbital.Orbital("NPP",
                              line1="1 37849U 11061A   12017.90990040 -.00000112  00000-0 -32693-4 0   772",
                              line2="2 37849  98.7026 317.8811 0001845  92.4533 267.6830 14.19582686 11574")
        dobj = datetime(2012, 1, 18, 8, 4, 19)
        orbnum = sat.get_orbit_number(dobj)
        self.assertEqual(orbnum, 1163)

    def test_sublonlat(self):
        sat = orbital.Orbital("ISS (ZARYA)",
                              line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
                              line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029")
        d = datetime(2003, 3, 23, 0, 3, 22)
        lon, lat, alt = sat.get_lonlatalt(d)
        expected_lon = -68.199894472013213
        expected_lat = 23.159747677881075
        expected_alt = 392.01953430856935
        self.failUnless(
            np.abs(lon - expected_lon) < eps_deg, 'Calculation of sublon failed')
        self.failUnless(
            np.abs(lat - expected_lat) < eps_deg, 'Calculation of sublat failed')
        self.failUnless(
            np.abs(alt - expected_alt) < eps_deg, 'Calculation of altitude failed')

    def test_observer_look(self):
        sat = orbital.Orbital("ISS (ZARYA)",
                              line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
                              line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029")
        d = datetime(2003, 3, 23, 0, 3, 22)
        az, el = sat.get_observer_look(d, -84.39733, 33.775867, 0)
        expected_az = 122.45169655331965
        expected_el = 1.9800219611255456
        self.failUnless(
            np.abs(az - expected_az) < eps_deg, 'Calculation of azimut failed')
        self.failUnless(
            np.abs(el - expected_el) < eps_deg, 'Calculation of elevation failed')


def suite():
    """The suite for orbital calculations"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestOrbital))

    return mysuite
