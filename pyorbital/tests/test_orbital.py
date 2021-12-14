#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012-2014 Martin Raspaud

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

"""Test the geoloc orbital.
"""

import unittest
try:
    from unittest import mock
except ImportError:
    import mock
from datetime import datetime, timedelta

import numpy as np

from pyorbital import orbital

eps_deg = 10e-3


class Test(unittest.TestCase):

    def test_get_orbit_number(self):
        """Testing getting the orbitnumber from the tle"""
        sat = orbital.Orbital("NPP",
                              line1="1 37849U 11061A   12017.90990040 "
                                    "-.00000112  00000-0 -32693-4 0   772",
                              line2="2 37849  98.7026 317.8811 0001845  "
                                    "92.4533 267.6830 14.19582686 11574")
        dobj = datetime(2012, 1, 18, 8, 4, 19)
        orbnum = sat.get_orbit_number(dobj)
        self.assertEqual(orbnum, 1163)

    def test_sublonlat(self):
        sat = orbital.Orbital("ISS (ZARYA)",
                              line1="1 25544U 98067A   03097.78853147  "
                                    ".00021906  00000-0  28403-3 0  8652",
                              line2="2 25544  51.6361  13.7980 0004256  "
                                    "35.6671  59.2566 15.58778559250029")
        d = datetime(2003, 3, 23, 0, 3, 22)
        lon, lat, alt = sat.get_lonlatalt(d)
        expected_lon = -68.199894472013213
        expected_lat = 23.159747677881075
        expected_alt = 392.01953430856935
        self.assertTrue(np.abs(lon - expected_lon) < eps_deg,
                        'Calculation of sublon failed')
        self.assertTrue(np.abs(lat - expected_lat) < eps_deg,
                        'Calculation of sublat failed')
        self.assertTrue(np.abs(alt - expected_alt) < eps_deg,
                        'Calculation of altitude failed')

    def test_observer_look(self):
        sat = orbital.Orbital("ISS (ZARYA)",
                              line1="1 25544U 98067A   03097.78853147  "
                                    ".00021906  00000-0  28403-3 0  8652",
                              line2="2 25544  51.6361  13.7980 0004256  "
                                    "35.6671  59.2566 15.58778559250029")
        d = datetime(2003, 3, 23, 0, 3, 22)
        az, el = sat.get_observer_look(d, -84.39733, 33.775867, 0)
        expected_az = 122.45169655331965
        expected_el = 1.9800219611255456
        self.assertTrue(np.abs(az - expected_az) < eps_deg,
                        'Calculation of azimut failed')
        self.assertTrue(np.abs(el - expected_el) < eps_deg,
                        'Calculation of elevation failed')

    def test_orbit_num_an(self):
        sat = orbital.Orbital("METOP-A",
                              line1="1 29499U 06044A   11254.96536486  "
                                    ".00000092  00000-0  62081-4 0  5221",
                              line2="2 29499  98.6804 312.6735 0001758 "
                                    "111.9178 248.2152 14.21501774254058")
        d = datetime(2011, 9, 14, 5, 30)
        self.assertEqual(sat.get_orbit_number(d), 25437)

    def test_orbit_num_non_an(self):
        sat = orbital.Orbital("METOP-A",
                              line1="1 29499U 06044A   13060.48822809  "
                                    ".00000017  00000-0  27793-4 0  9819",
                              line2="2 29499  98.6639 121.6164 0001449  "
                                    "71.9056  43.3132 14.21510544330271")
        dt = np.timedelta64(98, 'm')
        self.assertEqual(sat.get_orbit_number(sat.tle.epoch + dt), 33028)

    def test_orbit_num_equator(self):
        sat = orbital.Orbital("SUOMI NPP",
                              line1="1 37849U 11061A   13061.24611272  "
                                    ".00000048  00000-0  43679-4 0  4334",
                              line2="2 37849  98.7444   1.0588 0001264  "
                                    "63.8791 102.8546 14.19528338 69643")
        t1 = datetime(2013, 3, 2, 22, 2, 25)
        t2 = datetime(2013, 3, 2, 22, 2, 26)
        on1 = sat.get_orbit_number(t1)
        on2 = sat.get_orbit_number(t2)
        self.assertEqual(on1, 6973)
        self.assertEqual(on2, 6974)
        pos1, vel1 = sat.get_position(t1, normalize=False)
        pos2, vel2 = sat.get_position(t2, normalize=False)
        del vel1, vel2
        self.assertTrue(pos1[2] < 0)
        self.assertTrue(pos2[2] > 0)

    def test_get_next_passes_apogee(self):
        """Regression test #22."""
        line1 = "1 24793U 97020B   18065.48735489  " \
                ".00000075  00000-0  19863-4 0  9994"
        line2 = "2 24793  86.3994 209.3241 0002020  " \
                "89.8714 270.2713 14.34246429 90794"

        orb = orbital.Orbital('IRIDIUM 7 [+]', line1=line1, line2=line2)
        d = datetime(2018, 3, 7, 3, 30, 15)
        res = orb.get_next_passes(d, 1, 170.556, -43.368, 0.5, horizon=40)
        self.assertTrue(abs(
            res[0][2] - datetime(2018, 3, 7, 3, 48, 13, 178439)) <
            timedelta(seconds=0.01))

    def test_get_next_passes_tricky(self):
        """ Check issue #34 for reference """
        line1 = "1 43125U 18004Q   18251.42128650 " \
            "+.00001666 +00000-0 +73564-4 0  9991"

        line2 = "2 43125 097.5269 314.3317 0010735 "\
            "157.6344 202.5362 15.23132245036381"

        orb = orbital.Orbital('LEMUR-2-BROWNCOW', line1=line1, line2=line2)
        d = datetime(2018, 9, 8)

        res = orb.get_next_passes(d, 72, -8.174163, 51.953319, 0.05, horizon=5)

        self.assertTrue(abs(
            res[0][2] - datetime(2018, 9, 8, 9, 5, 46, 375248)) <
            timedelta(seconds=0.01))
        self.assertTrue(abs(
            res[-1][2] - datetime(2018, 9, 10, 22, 15, 3, 143469)) <
            timedelta(seconds=0.01))

        self.assertTrue(len(res) == 15)

    def test_get_next_passes_issue_22(self):
        """Check that max"""
        line1 = '1 28654U 05018A   21083.16603416  .00000102  00000-0  79268-4 0  9999'
        line2 = '2 28654  99.0035 147.6583 0014816 159.4931 200.6838 14.12591533816498'

        orb = orbital.Orbital("NOAA 18", line1=line1, line2=line2)
        t = datetime(2021, 3, 9, 22)
        next_passes = orb.get_next_passes(t, 1, -15.6335, 27.762, 0.)
        rise, fall, max_elevation = next_passes[0]
        assert rise < max_elevation < fall
        print(next_passes)

    @mock.patch('pyorbital.orbital.Orbital.get_lonlatalt')
    def test_utc2local(self, get_lonlatalt):
        get_lonlatalt.return_value = -45, None, None
        sat = orbital.Orbital("METOP-A",
                              line1="1 29499U 06044A   13060.48822809  "
                                    ".00000017  00000-0  27793-4 0  9819",
                              line2="2 29499  98.6639 121.6164 0001449  "
                                    "71.9056  43.3132 14.21510544330271")
        self.assertEqual(sat.utc2local(datetime(2009, 7, 1, 12)),
                         datetime(2009, 7, 1, 9))

    @mock.patch('pyorbital.orbital.Orbital.utc2local')
    @mock.patch('pyorbital.orbital.Orbital.get_orbit_number')
    def test_get_equatorial_crossing_time(self, get_orbit_number, utc2local):
        def get_orbit_number_patched(utc_time, **kwargs):
            utc_time = np.datetime64(utc_time)
            diff = (utc_time - np.datetime64('2009-07-01 12:38:12')) / np.timedelta64(7200, 's')
            return 1234 + diff

        get_orbit_number.side_effect = get_orbit_number_patched
        utc2local.return_value = 'local_time'
        sat = orbital.Orbital("METOP-A",
                              line1="1 29499U 06044A   13060.48822809  "
                                    ".00000017  00000-0  27793-4 0  9819",
                              line2="2 29499  98.6639 121.6164 0001449  "
                                    "71.9056  43.3132 14.21510544330271")

        # Ascending node
        res = sat.get_equatorial_crossing_time(tstart=datetime(2009, 7, 1, 12),
                                               tend=datetime(2009, 7, 1, 13))
        exp = datetime(2009, 7, 1, 12, 38, 12)
        self.assertTrue((res - exp) < timedelta(seconds=0.01))

        # Descending node
        res = sat.get_equatorial_crossing_time(tstart=datetime(2009, 7, 1, 12),
                                               tend=datetime(2009, 7, 1, 14, 0),
                                               node='descending')
        exp = datetime(2009, 7, 1, 13, 38, 12)
        self.assertTrue((res - exp) < timedelta(seconds=0.01))

        # Conversion to local time
        res = sat.get_equatorial_crossing_time(tstart=datetime(2009, 7, 1, 12),
                                               tend=datetime(2009, 7, 1, 14),
                                               local_time=True)
        self.assertEqual(res, 'local_time')


class TestGetObserverLook(unittest.TestCase):
    """Test the get_observer_look function"""

    def setUp(self):
        self.t = datetime(2018, 1, 1, 0, 0, 0)
        self.sat_lon = np.array([[-89.5, -89.4, -89.5, -89.4],
                                 [-89.3, -89.2, -89.3, -89.2]])
        self.sat_lat = np.array([[45.5, 45.4, 45.5, 45.4],
                                 [45.3, 40.2, 45.3, 40.2]])
        self.sat_alt = 35786 * np.ones((2, 4))
        self.lon = np.array([[-85.5, -85.4, -89.5, -99.4],
                             [-85.3, -89.2, -89.3, -79.2]])
        self.lat = np.array([[40.5, 40.4, 65.5, 45.4],
                             [40.3, 40.2, 25.3, 40.2]])
        self.alt = np.zeros((2, 4))
        self.exp_azi = np.array([[331.00275902, 330.95954165, 180, 86.435411],
                                 [330.91642994, 180, 0, 273.232073]])
        self.exp_elev = np.array([[83.18070976, 83.17788976, 66.548467, 81.735221],
                                  [83.17507167, 90, 66.559906, 81.010018]])

    def test_basic_numpy(self):
        """Test with numpy array inputs"""
        from pyorbital import orbital
        azi, elev = orbital.get_observer_look(self.sat_lon, self.sat_lat,
                                              self.sat_alt, self.t,
                                              self.lon, self.lat, self.alt)
        np.testing.assert_allclose(azi, self.exp_azi)
        np.testing.assert_allclose(elev, self.exp_elev)

    def test_basic_dask(self):
        """Test with dask array inputs"""
        from pyorbital import orbital
        import dask.array as da
        sat_lon = da.from_array(self.sat_lon, chunks=2)
        sat_lat = da.from_array(self.sat_lat, chunks=2)
        sat_alt = da.from_array(self.sat_alt, chunks=2)
        lon = da.from_array(self.lon, chunks=2)
        lat = da.from_array(self.lat, chunks=2)
        alt = da.from_array(self.alt, chunks=2)
        azi, elev = orbital.get_observer_look(sat_lon, sat_lat,
                                              sat_alt, self.t,
                                              lon, lat, alt)
        np.testing.assert_allclose(azi.compute(), self.exp_azi)
        np.testing.assert_allclose(elev.compute(), self.exp_elev)

    def test_xarray_with_numpy(self):
        """Test with xarray DataArray with numpy array as inputs"""
        from pyorbital import orbital
        import xarray as xr

        def _xarr_conv(input):
            return xr.DataArray(input)
        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = orbital.get_observer_look(sat_lon, sat_lat,
                                              sat_alt, self.t,
                                              lon, lat, alt)
        np.testing.assert_allclose(azi.data, self.exp_azi)
        np.testing.assert_allclose(elev.data, self.exp_elev)

    def test_xarray_with_dask(self):
        """Test with xarray DataArray with dask array as inputs"""
        from pyorbital import orbital
        import dask.array as da
        import xarray as xr

        def _xarr_conv(input):
            return xr.DataArray(da.from_array(input, chunks=2))
        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = orbital.get_observer_look(sat_lon, sat_lat,
                                              sat_alt, self.t,
                                              lon, lat, alt)
        np.testing.assert_allclose(azi.data.compute(), self.exp_azi)
        np.testing.assert_allclose(elev.data.compute(), self.exp_elev)

    def test_scalar(self):
        """Test with scalar inputs."""
        from pyorbital.orbital import get_observer_look
        (azi, elev) = get_observer_look(0, 0, 30_000_000, self.t, 0, 0, 0)
        np.testing.assert_allclose(elev, 90)


class TestGetObserverLookNadir(unittest.TestCase):
    """Test the get_observer_look function when satellite is at nadir."""

    def setUp(self):
        """Setup for test observer at nadir.
        Note that rounding error differs between array types.
        With 1000 elements a test gives:
        1 error for basic numpy
        41 errors for basic dask
        63 errors for xarray with dask
        2 error for xarray with numpy
        """
        rng = np.random.RandomState(125)
        self.t = datetime(2018, 1, 1, 0, 0, 0)
        self.sat_lon = 360 * rng.rand(100) - 180
        self.sat_lat = 180 * rng.rand(100) - 90
        self.sat_alt = rng.rand(100) + 850
        self.lon = self.sat_lon  # + 10E-17
        self.lat = self.sat_lat  # + 10E-17
        self.alt = np.zeros((100))
        self.exp_elev = np.zeros((100)) + 90

    def test_basic_numpy(self):
        """Test with numpy array inputs"""
        from pyorbital import orbital
        azi, elev = orbital.get_observer_look(self.sat_lon, self.sat_lat,
                                              self.sat_alt, self.t,
                                              self.lon, self.lat, self.alt)
        self.assertEqual(np.sum(np.isnan(azi)), 0)
        self.assertFalse(np.isnan(azi).any())
        np.testing.assert_allclose(elev, self.exp_elev)

    def test_basic_dask(self):
        """Test with dask array inputs"""
        from pyorbital import orbital
        import dask.array as da
        sat_lon = da.from_array(self.sat_lon, chunks=2)
        sat_lat = da.from_array(self.sat_lat, chunks=2)
        sat_alt = da.from_array(self.sat_alt, chunks=2)
        lon = da.from_array(self.lon, chunks=2)
        lat = da.from_array(self.lat, chunks=2)
        alt = da.from_array(self.alt, chunks=2)
        azi, elev = orbital.get_observer_look(sat_lon, sat_lat,
                                              sat_alt, self.t,
                                              lon, lat, alt)
        self.assertEqual(np.sum(np.isnan(azi)), 0)
        self.assertFalse(np.isnan(azi).any())
        np.testing.assert_allclose(elev.compute(), self.exp_elev)

    def test_xarray_with_numpy(self):
        """Test with xarray DataArray with numpy array as inputs"""
        from pyorbital import orbital
        import xarray as xr

        def _xarr_conv(input):
            return xr.DataArray(input)
        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = orbital.get_observer_look(sat_lon, sat_lat,
                                              sat_alt, self.t,
                                              lon, lat, alt)
        self.assertEqual(np.sum(np.isnan(azi)), 0)
        self.assertFalse(np.isnan(azi).any())
        np.testing.assert_allclose(elev.data, self.exp_elev)

    def test_xarray_with_dask(self):
        """Test with xarray DataArray with dask array as inputs"""
        from pyorbital import orbital
        import dask.array as da
        import xarray as xr

        def _xarr_conv(input):
            return xr.DataArray(da.from_array(input, chunks=2))
        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = orbital.get_observer_look(sat_lon, sat_lat,
                                              sat_alt, self.t,
                                              lon, lat, alt)
        self.assertEqual(np.sum(np.isnan(azi)), 0)
        self.assertFalse(np.isnan(azi).any())
        np.testing.assert_allclose(elev.data.compute(), self.exp_elev)


class TestRegressions(unittest.TestCase):
    """Test regressions."""

    def test_63(self):
        """Check that no runtimewarning is raised, #63."""
        import warnings
        from pyorbital.orbital import Orbital
        from dateutil import parser
        warnings.filterwarnings('error')
        orb = Orbital("Suomi-NPP",
                      line1="1 37849U 11061A   19292.84582509  .00000011  00000-0  25668-4 0  9997",
                      line2="2 37849  98.7092 229.3263 0000715  98.5313 290.6262 14.19554485413345")
        orb.get_next_passes(parser.parse("2019-10-21 16:00:00"), 12, 123.29736, -13.93763, 0)
        warnings.filterwarnings('default')


def suite():
    """The suite for test_orbital
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(Test))
    mysuite.addTest(loader.loadTestsFromTestCase(TestGetObserverLook))
    mysuite.addTest(loader.loadTestsFromTestCase(TestGetObserverLookNadir))
    mysuite.addTest(loader.loadTestsFromTestCase(TestRegressions))

    return mysuite
