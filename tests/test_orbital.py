import unittest
from datetime import datetime, timedelta

import numpy as np

from pyorbital import orbital

eps_deg = 10e-3

class Test(unittest.TestCase):

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
        self.failUnless(np.abs(lon - expected_lon) < eps_deg, 'Calculation of sublon failed')
        self.failUnless(np.abs(lat - expected_lat) < eps_deg, 'Calculation of sublat failed')
        self.failUnless(np.abs(alt - expected_alt) < eps_deg, 'Calculation of altitude failed')
        

    def test_observer_look(self):
        sat = orbital.Orbital("ISS (ZARYA)", 
            line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652", 
            line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029")
        d = datetime(2003, 3, 23, 0, 3, 22)
        az, el = sat.get_observer_look(d, -84.39733, 33.775867, 0)
        expected_az = 122.45169655331965
        expected_el = 1.9800219611255456
        self.failUnless(np.abs(az - expected_az) < eps_deg, 'Calculation of azimut failed')
        self.failUnless(np.abs(el - expected_el) < eps_deg, 'Calculation of elevation failed')
        
    def test_orbit_num_an(self):
        sat = orbital.Orbital("METOP-A", 
            line1="1 29499U 06044A   11254.96536486  .00000092  00000-0  62081-4 0  5221", 
            line2="2 29499  98.6804 312.6735 0001758 111.9178 248.2152 14.21501774254058")
        d = datetime(2011, 9, 14, 5, 30)
        self.assertEqual(sat.get_orbit_number(d), 25437)
        
    def test_orbit_num_non_an(self):
        sat = orbital.Orbital("METOP-A", 
            line1="1 29499U 06044A   13060.48822809  .00000017  00000-0  27793-4 0  9819", 
            line2="2 29499  98.6639 121.6164 0001449  71.9056  43.3132 14.21510544330271")
        dt = timedelta(minutes=98)
        self.assertEqual(sat.get_orbit_number(sat.tle.epoch + dt), 33028)
        
    def test_orbit_num_equator(self):
        sat = orbital.Orbital("SUOMI NPP", 
            line1="1 37849U 11061A   13061.24611272  .00000048  00000-0  43679-4 0  4334", 
            line2="2 37849  98.7444   1.0588 0001264  63.8791 102.8546 14.19528338 69643")
        t1 = datetime(2013, 3, 2, 22, 2, 25)
        t2 = datetime(2013, 3, 2, 22, 2, 26)
        on1 = sat.get_orbit_number(t1)
        on2 = sat.get_orbit_number(t2)
        self.assertEqual(on1, 6973)
        self.assertEqual(on2, 6974)
        pos1, vel1 = sat.get_position(t1, normalize=False)
        pos2, vel2 = sat.get_position(t2, normalize=False)
        self.assertTrue(pos1[2] < 0)
        self.assertTrue(pos2[2] > 0)
        
if __name__ == "__main__":
    unittest.main()
