#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Martin Raspaud
#
# Author(s):
#
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Test TLE file reading, TLE downloading and stroging TLEs to database."""


from pyorbital.tlefile import Tle
import datetime
import unittest
from unittest import mock

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
        """Check the *tle* instance against predetermined values."""
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
        """Test parsing from line elements."""
        tle = Tle("ISS (ZARYA)", line1=line1, line2=line2)
        self.check_example(tle)

    def test_from_file(self):
        """Test reading and parsing from a file."""
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


class TestDownloader(unittest.TestCase):
    """Test TLE downloader."""

    def setUp(self):
        """Create a downloader instance."""
        from pyorbital.tlefile import Downloader
        self.config = {}
        self.dl = Downloader(self.config)

    def test_init(self):
        """Test the initialization."""
        assert self.dl.config is self.config

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_plain_tle(self, requests):
        """Test downloading and a TLE file from internet."""
        requests.get = mock.MagicMock()
        # The return value of requests.get()
        req = mock.MagicMock()
        req.status_code = 200
        req.text = '\n'.join((line0, line1, line2))
        requests.get.return_value = req

        # Not configured
        self.dl.config["downloaders"] = {}
        res = self.dl.fetch_plain_tle()
        self.assertTrue(res == {})
        requests.get.assert_not_called()

        # Two sources, one with multiple locations
        self.dl.config["downloaders"] = {
            "fetch_plain_tle": {
                "source_1": ["mocked_url_1", "mocked_url_2", "mocked_url_3"],
                "source_2": ["mocked_url_4"]
            }
        }
        res = self.dl.fetch_plain_tle()
        self.assertTrue("source_1" in res)
        self.assertEqual(len(res["source_1"]), 3)
        self.assertEqual(res["source_1"][0].line1, line1)
        self.assertEqual(res["source_1"][0].line2, line2)
        self.assertTrue("source_2" in res)
        self.assertEqual(len(res["source_2"]), 1)
        self.assertTrue(mock.call("mocked_url_1") in requests.get.mock_calls)
        self.assertEqual(len(requests.get.mock_calls), 4)

        # Reset mocks
        requests.get.reset_mock()
        req.reset_mock()

        # No data returned because the server is a teapot
        req.status_code = 418
        res = self.dl.fetch_plain_tle()
        # The sources are in the dict ...
        self.assertEqual(len(res), 2)
        # ... but there are no TLEs
        self.assertEqual(len(res["source_1"]), 0)
        self.assertEqual(len(res["source_2"]), 0)
        self.assertTrue(mock.call("mocked_url_1") in requests.get.mock_calls)
        self.assertEqual(len(requests.get.mock_calls), 4)

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_spacetrack(self, requests):
        """Test downloading and TLEs from space-track.org."""
        mock_post = mock.MagicMock()
        mock_get = mock.MagicMock()
        mock_session = mock.MagicMock()
        mock_session.post = mock_post
        mock_session.get = mock_get
        requests.Session.return_value.__enter__.return_value = mock_session

        tle_text = '\n'.join((line0, line1, line2))
        self.dl.config["platforms"] = {
            25544: 'ISS'
        }
        self.dl.config["downloaders"] = {
            "fetch_spacetrack": {
                "user": "username",
                "password": "passw0rd"
            }
        }

        # Login fails, because the server is a teapot
        mock_post.return_value.status_code = 418
        res = self.dl.fetch_spacetrack()
        # Empty list of TLEs is returned
        self.assertTrue(res == [])
        # The login was anyway attempted
        mock_post.assert_called_with(
            'https://www.space-track.org/ajaxauth/login',
            data={'identity': 'username', 'password': 'passw0rd'})

        # Login works, but something is wrong (teapot) when asking for data
        mock_post.return_value.status_code = 200
        mock_get.return_value.status_code = 418
        res = self.dl.fetch_spacetrack()
        self.assertTrue(res == [])
        mock_get.assert_called_with("https://www.space-track.org/"
                                    "basicspacedata/query/class/tle_latest/"
                                    "ORDINAL/1/NORAD_CAT_ID/25544/format/tle")

        # Data is received
        mock_get.return_value.status_code = 200
        mock_get.return_value.text = tle_text
        res = self.dl.fetch_spacetrack()
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0].line1, line1)
        self.assertEqual(res[0].line2, line2)

    def test_read_tle_files(self):
        """Test reading TLE files from a file system."""
        from tempfile import TemporaryDirectory
        import os

        tle_text = '\n'.join((line0, line1, line2))

        save_dir = TemporaryDirectory()
        with save_dir:
            fname = os.path.join(save_dir.name, 'tle_20200129_1600.txt')
            with open(fname, 'w') as fid:
                fid.write(tle_text)
            # Add a non-existent file, it shouldn't cause a crash
            nonexistent = os.path.join(save_dir.name, 'not_here.txt')
            # Use a wildcard to collect files (passed to glob)
            starred_fname = os.path.join(save_dir.name, 'tle*txt')
            self.dl.config["downloaders"] = {
                "read_tle_files": {
                    "paths": [fname, nonexistent, starred_fname]
                }
            }
            res = self.dl.read_tle_files()
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0].line1, line1)
        self.assertEqual(res[0].line2, line2)


def suite():
    """Create the test suite for test_tlefile."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TLETest))
    mysuite.addTest(loader.loadTestsFromTestCase(TestDownloader))

    return mysuite
