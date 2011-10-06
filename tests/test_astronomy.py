import unittest

from datetime import datetime
import numpy as np
import pyorbital.astronomy as astr



class TestAstronomy(unittest.TestCase):

    def setUp(self):
        pass

    def test_jdays(self):
        """Test julian day functions.
        """

        t = datetime(2000, 1, 1, 12, 0)
        self.assertEqual(astr.jdays(t), 2451545.0)
        self.assertEqual(astr.jdays2000(t), 0)
        t = datetime(2009, 10, 8, 14, 30)
        self.assertEqual(astr.jdays(t), 2455113.1041666665)
        self.assertEqual(astr.jdays2000(t), 3568.1041666666665)
        
    def test_sunangles(self):
        """Test the sun-angle calculations:
        """
        lat, lon = 58.6167, 16.1833 # Norrkoping
        time_slot = datetime(2011, 9, 23, 12, 0)
        
        sun_theta = astr.sun_zenith_angle(time_slot, lon, lat)
        self.assertEqual(sun_theta, 60.371433482557833)
        sun_theta = astr.sun_zenith_angle(time_slot, 0., 0.)
        self.assertEqual(sun_theta, 1.8751916863323426)

        
if __name__ == '__main__':
    unittest.main()

