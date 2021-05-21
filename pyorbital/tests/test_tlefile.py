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
import os
from contextlib import suppress

line0 = "ISS (ZARYA)"
line1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"

line1_2 = "1 38771U 12049A   21137.30264622  .00000000  00000+0 -49996-5 0 00017"
line2_2 = "2 38771  98.7162 197.7716 0002383 106.1049 122.6344 14.21477797449453"

tle_xml = '\n'.join(
    ('<?xml version="1.0" encoding="UTF-8"?>',
        '<multi-mission-administrative-message>',
        '<message>',
        '<two-line-elements>',
        '<navigation>',
        '<line-1>' + line1 + '</line-1>',
        '<line-2>' + line2 + '</line-2>',
        '</navigation>',
        '</two-line-elements>',
        '</message>',
        '<message>',
        '<two-line-elements>',
        '<navigation>',
        '<line-1>' + line1_2 + '</line-1>',
        '<line-2>' + line2_2 + '</line-2>',
        '</navigation>',
        '</two-line-elements>',
        '</message>',
        '</multi-mission-administrative-message>'))


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

    def test_from_mmam_xml(self):
        """Test reading from an MMAM XML file."""
        from tempfile import TemporaryDirectory

        save_dir = TemporaryDirectory()
        with save_dir:
            fname = os.path.join(save_dir.name, '20210420_Metop-B_ADMIN_MESSAGE_NO_127.xml')
            with open(fname, 'w') as fid:
                fid.write(tle_xml)
            tle = Tle("", tle_file=fname)
        self.check_example(tle)


FETCH_PLAIN_TLE_CONFIG = {
    "fetch_plain_tle": {
        "source_1": ["mocked_url_1", "mocked_url_2", "mocked_url_3"],
        "source_2": ["mocked_url_4"]
    }
}
FETCH_SPACETRACK_CONFIG = {
    "fetch_spacetrack": {
        "user": "username",
        "password": "passw0rd"
    }
}


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
    def test_fetch_plain_tle_not_configured(self, requests):
        """Test downloading and a TLE file from internet."""
        requests.get = mock.MagicMock()
        requests.get.return_value = _get_req_response(200)

        # Not configured
        self.dl.config["downloaders"] = {}
        res = self.dl.fetch_plain_tle()
        self.assertTrue(res == {})
        requests.get.assert_not_called()

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_plain_tle_two_sources(self, requests):
        """Test downloading and a TLE file from internet."""
        requests.get = mock.MagicMock()
        requests.get.return_value = _get_req_response(200)

        # Two sources, one with multiple locations
        self.dl.config["downloaders"] = FETCH_PLAIN_TLE_CONFIG

        res = self.dl.fetch_plain_tle()
        self.assertTrue("source_1" in res)
        self.assertEqual(len(res["source_1"]), 3)
        self.assertEqual(res["source_1"][0].line1, line1)
        self.assertEqual(res["source_1"][0].line2, line2)
        self.assertTrue("source_2" in res)
        self.assertEqual(len(res["source_2"]), 1)
        self.assertTrue(mock.call("mocked_url_1") in requests.get.mock_calls)
        self.assertEqual(len(requests.get.mock_calls), 4)

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_plain_tle_server_is_a_teapot(self, requests):
        """Test downloading and a TLE file from internet."""
        requests.get = mock.MagicMock()
        # No data returned because the server is a teapot
        requests.get.return_value = _get_req_response(418)

        # Two sources, one with multiple locations
        self.dl.config["downloaders"] = FETCH_PLAIN_TLE_CONFIG

        res = self.dl.fetch_plain_tle()
        # The sources are in the dict ...
        self.assertEqual(len(res), 2)
        # ... but there are no TLEs
        self.assertEqual(len(res["source_1"]), 0)
        self.assertEqual(len(res["source_2"]), 0)
        self.assertTrue(mock.call("mocked_url_1") in requests.get.mock_calls)
        self.assertEqual(len(requests.get.mock_calls), 4)

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_spacetrack_login_fails(self, requests):
        """Test downloading and TLEs from space-track.org."""
        mock_post = mock.MagicMock()
        mock_session = mock.MagicMock()
        mock_session.post = mock_post
        requests.Session.return_value.__enter__.return_value = mock_session

        self.dl.config["platforms"] = {
            25544: 'ISS'
        }
        self.dl.config["downloaders"] = FETCH_SPACETRACK_CONFIG

        # Login fails, because the server is a teapot
        mock_post.return_value.status_code = 418
        res = self.dl.fetch_spacetrack()
        # Empty list of TLEs is returned
        self.assertTrue(res == [])
        # The login was anyway attempted
        mock_post.assert_called_with(
            'https://www.space-track.org/ajaxauth/login',
            data={'identity': 'username', 'password': 'passw0rd'})

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_spacetrack_get_fails(self, requests):
        """Test downloading and TLEs from space-track.org."""
        mock_post = mock.MagicMock()
        mock_get = mock.MagicMock()
        mock_session = mock.MagicMock()
        mock_session.post = mock_post
        mock_session.get = mock_get
        requests.Session.return_value.__enter__.return_value = mock_session

        self.dl.config["platforms"] = {
            25544: 'ISS'
        }
        self.dl.config["downloaders"] = FETCH_SPACETRACK_CONFIG

        # Login works, but something is wrong (teapot) when asking for data
        mock_post.return_value.status_code = 200
        mock_get.return_value.status_code = 418
        res = self.dl.fetch_spacetrack()
        self.assertTrue(res == [])
        mock_get.assert_called_with("https://www.space-track.org/"
                                    "basicspacedata/query/class/tle_latest/"
                                    "ORDINAL/1/NORAD_CAT_ID/25544/format/tle")

    @mock.patch('pyorbital.tlefile.requests')
    def test_fetch_spacetrack_success(self, requests):
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
        self.dl.config["downloaders"] = FETCH_SPACETRACK_CONFIG

        # Login works and data is received
        mock_post.return_value.status_code = 200
        mock_get.return_value.status_code = 200
        mock_get.return_value.text = tle_text
        res = self.dl.fetch_spacetrack()
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0].line1, line1)
        self.assertEqual(res[0].line2, line2)

    def test_read_tle_files(self):
        """Test reading TLE files from a file system."""
        from tempfile import TemporaryDirectory

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

    def test_read_xml_admin_messages(self):
        """Test reading TLE files from a file system."""
        from tempfile import TemporaryDirectory

        save_dir = TemporaryDirectory()
        with save_dir:
            fname = os.path.join(save_dir.name, '20210420_Metop-B_ADMIN_MESSAGE_NO_127.xml')
            with open(fname, 'w') as fid:
                fid.write(tle_xml)
            # Add a non-existent file, it shouldn't cause a crash
            nonexistent = os.path.join(save_dir.name, 'not_here.txt')
            # Use a wildcard to collect files (passed to glob)
            starred_fname = os.path.join(save_dir.name, '*.xml')
            self.dl.config["downloaders"] = {
                "read_xml_admin_messages": {
                    "paths": [fname, nonexistent, starred_fname]
                }
            }
            res = self.dl.read_xml_admin_messages()

        # There are two sets of TLEs in the file.  And as the same file is
        # parsed twice, 4 TLE objects are returned
        self.assertEqual(len(res), 4)
        self.assertEqual(res[0].line1, line1)
        self.assertEqual(res[0].line2, line2)
        self.assertEqual(res[1].line1, line1_2)
        self.assertEqual(res[1].line2, line2_2)


def _get_req_response(code):
    req = mock.MagicMock()
    req.status_code = code
    req.text = '\n'.join((line0, line1, line2))
    return req


class TestSQLiteTLE(unittest.TestCase):
    """Test saving TLE data to a SQLite database."""

    def setUp(self):
        """Create a database instance."""
        from pyorbital.tlefile import SQLiteTLE
        from pyorbital.tlefile import Tle
        from tempfile import TemporaryDirectory

        self.temp_dir = TemporaryDirectory()
        self.db_fname = os.path.join(self.temp_dir.name, 'tle.db')
        self.platforms = {25544: "ISS"}
        self.writer_config = {
            "output_dir": os.path.join(self.temp_dir.name, 'tle_dir'),
            "filename_pattern": "tle_%Y%m%d_%H%M%S.%f.txt",
            "write_name": True,
            "write_always": False
        }
        self.db = SQLiteTLE(self.db_fname, self.platforms, self.writer_config)
        self.tle = Tle('ISS', line1=line1, line2=line2)

    def tearDown(self):
        """Clean temporary files."""
        with suppress(PermissionError, NotADirectoryError):
            self.temp_dir.cleanup()

    def test_init(self):
        """Test that the init did what it should have."""
        from pyorbital.tlefile import table_exists, PLATFORM_NAMES_TABLE

        columns = [col.strip() for col in
                   PLATFORM_NAMES_TABLE.strip('()').split(',')]
        num_columns = len(columns)

        self.assertTrue(os.path.exists(self.db_fname))
        self.assertTrue(table_exists(self.db.db, "platform_names"))
        res = self.db.db.execute('select * from platform_names')
        names = [description[0] for description in res.description]
        self.assertEqual(len(names), num_columns)
        for col in columns:
            self.assertTrue(col.split(' ')[0] in names)

    def test_update_db(self):
        """Test updating database with new data."""
        from pyorbital.tlefile import (table_exists, SATID_TABLE,
                                       ISO_TIME_FORMAT)

        # Get the column names
        columns = [col.strip() for col in
                   SATID_TABLE.replace("'{}' (", "").strip(')').split(',')]
        # Platform number
        satid = str(list(self.platforms.keys())[0])

        # Data from a platform that isn't configured
        self.db.platforms = {}
        self.db.update_db(self.tle, 'foo')
        self.assertFalse(table_exists(self.db.db, satid))
        self.assertFalse(self.db.updated)

        # Configured platform
        self.db.platforms = self.platforms
        self.db.update_db(self.tle, 'foo')
        self.assertTrue(table_exists(self.db.db, satid))
        self.assertTrue(self.db.updated)

        # Check that all the columns were added
        res = self.db.db.execute("select * from '%s'" % satid)
        names = [description[0] for description in res.description]
        for col in columns:
            self.assertTrue(col.split(' ')[0] in names)

        # Check the data
        data = res.fetchall()
        self.assertEqual(len(data), 1)
        # epoch
        self.assertEqual(data[0][0], '2008-09-20T12:25:40.104192')
        # TLE
        self.assertEqual(data[0][1], '\n'.join((line1, line2)))
        # Date when the data were added should be close to current time
        date_added = datetime.datetime.strptime(data[0][2], ISO_TIME_FORMAT)
        now = datetime.datetime.utcnow()
        self.assertTrue((now - date_added).total_seconds() < 1.0)
        # Source of the data
        self.assertTrue(data[0][3] == 'foo')

        # Try to add the same data again. Nothing should change even
        # if the source is different if the epoch is the same
        self.db.update_db(self.tle, 'bar')
        res = self.db.db.execute("select * from '%s'" % satid)
        data = res.fetchall()
        self.assertEqual(len(data), 1)
        date_added2 = datetime.datetime.strptime(data[0][2], ISO_TIME_FORMAT)
        self.assertEqual(date_added, date_added2)
        # Source of the data
        self.assertTrue(data[0][3] == 'foo')

    def test_write_tle_txt(self):
        """Test reading data from the database and writing it to a file."""
        import glob
        tle_dir = self.writer_config["output_dir"]

        # Put some data in the database
        self.db.update_db(self.tle, 'foo')

        # Fake that the database hasn't been updated
        self.db.updated = False

        # Try to dump the data to disk
        self.db.write_tle_txt()

        # The output dir hasn't been created
        self.assertFalse(os.path.exists(tle_dir))

        self.db.updated = True
        self.db.write_tle_txt()

        # The dir should be there
        self.assertTrue(os.path.exists(tle_dir))
        # There should be one file in the directory
        files = glob.glob(os.path.join(tle_dir, 'tle_*txt'))
        self.assertEqual(len(files), 1)
        # The file should have been named with the date ('%' characters
        # not there anymore)
        self.assertTrue('%' not in files[0])
        # The satellite name should be in the file
        with open(files[0], 'r') as fid:
            data = fid.read().split('\n')
        self.assertEqual(len(data), 3)
        self.assertTrue('ISS' in data[0])
        self.assertEqual(data[1], line1)
        self.assertEqual(data[2], line2)

        # Call the writing again, nothing should be written. In
        # real-life this assumes a re-run has been done without new
        # TLE data
        self.db.updated = False
        self.db.write_tle_txt()
        files = glob.glob(os.path.join(tle_dir, 'tle_*txt'))
        self.assertEqual(len(files), 1)

        # Force writing with every call
        # Do not write the satellite name
        self.db.writer_config["write_always"] = True
        self.db.writer_config["write_name"] = False
        self.db.write_tle_txt()
        files = sorted(glob.glob(os.path.join(tle_dir, 'tle_*txt')))
        self.assertEqual(len(files), 2)
        with open(files[1], 'r') as fid:
            data = fid.read().split('\n')
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0], line1)
        self.assertEqual(data[1], line2)


def suite():
    """Create the test suite for test_tlefile."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TLETest))
    mysuite.addTest(loader.loadTestsFromTestCase(TestDownloader))
    mysuite.addTest(loader.loadTestsFromTestCase(TestSQLiteTLE))

    return mysuite
