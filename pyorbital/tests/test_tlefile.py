#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014-2024 Pytroll Community
#
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


import datetime
import logging
import os
import time
import unittest
from contextlib import suppress
from pathlib import Path
from tempfile import mkstemp
from unittest import mock

import pytest

LINE0 = "ISS (ZARYA)"
LINE1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
LINE2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"
LINE1_2 = "1 38771U 12049A   21137.30264622  .00000000  00000+0 -49996-5 0 00017"
LINE2_2 = "2 38771  98.7162 197.7716 0002383 106.1049 122.6344 14.21477797449453"

def _write_fake_platforms_txt_file(platforms_filename) -> None:
    with open(platforms_filename, "w") as platforms_file:
        platforms_file.write("""NOAA-18 28654
NOAA-19 33591
NOAA-20 43013
NOAA-21 54234
# ISS 25544
""")


# NOAA 18
# 1 28654U 05018A   23045.48509621  .00000446  00000+0  26330-3 0  9998
# 2 28654  98.9223 120.4228 0014233  11.3574 348.7916 14.12862494914152

def _write_fake_tle_file(tlefilename: Path) -> None:
    with open(tlefilename, "w") as tle_file:
        tle_file.write("""NOAA 20
1 43013U 17073A   23045.54907786  .00000253  00000+0  14081-3 0  9995
2 43013  98.7419 345.5839 0001610  80.3742 279.7616 14.19558274271576
NOAA 21 (JPSS-2)
1 54234U 22150A   23045.56664999  .00000332  00000+0  17829-3 0  9993
2 54234  98.7059 345.5113 0001226  81.6523 278.4792 14.19543871 13653
ISS (ZARYA)
1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
""")


@pytest.fixture
def fake_platforms_txt_file(tmp_path: Path) -> Path:
    """Make fake platforms.txt file."""
    filename = tmp_path / "platforms.txt"
    _write_fake_platforms_txt_file(filename)
    return filename


@pytest.fixture
def fake_tlefile(tmp_path: Path) -> Path:
    """Make fake tle file."""
    filename = tmp_path / "sometlefile.txt"
    _write_fake_tle_file(filename)
    return filename




NOAA19_2LINES = """1 33591U 09005A   21355.91138073  .00000074  00000+0  65091-4 0  9998
2 33591  99.1688  21.1338 0013414 329.8936  30.1462 14.12516400663123
"""
NOAA19_3LINES = "NOAA 19\n" + NOAA19_2LINES


tle_xml = "\n".join(
    ('<?xml version="1.0" encoding="UTF-8"?>',
        "<multi-mission-administrative-message>",
        "<message>",
        "<two-line-elements>",
        "<navigation>",
        "<line-1>" + LINE1 + "</line-1>",
        "<line-2>" + LINE2 + "</line-2>",
        "</navigation>",
        "</two-line-elements>",
        "</message>",
        "<message>",
        "<two-line-elements>",
        "<navigation>",
        "<line-1>" + LINE1_2 + "</line-1>",
        "<line-2>" + LINE2_2 + "</line-2>",
        "</navigation>",
        "</two-line-elements>",
        "</message>",
        "</multi-mission-administrative-message>"))


def test_read_tlefile_standard_platform_name(monkeypatch, fake_platforms_txt_file, fake_tlefile):
    """Test create a tle-object by reading tle data from file.

    Use Oscar naming matching name in platforms.txt.
    """
    from pyorbital import tlefile

    path_to_platforms_txt_file = fake_platforms_txt_file.parent
    monkeypatch.setenv("PYORBITAL_CONFIG_PATH", str(path_to_platforms_txt_file))

    tle_n21 = tlefile.read("NOAA-21", str(fake_tlefile))
    assert tle_n21.line1 == "1 54234U 22150A   23045.56664999  .00000332  00000+0  17829-3 0  9993"
    assert tle_n21.line2 == "2 54234  98.7059 345.5113 0001226  81.6523 278.4792 14.19543871 13653"


def test_read_tlefile_non_standard_platform_name(monkeypatch, fake_platforms_txt_file, fake_tlefile):
    """Test create a tle-object by reading tle data from file.

    Use naming matching what is in the TLE files, but non-standard (non Oscar) naming.
    """
    from pyorbital import tlefile

    path_to_platforms_txt_file = fake_platforms_txt_file.parent
    monkeypatch.setenv("PYORBITAL_CONFIG_PATH", str(path_to_platforms_txt_file))

    tle_n20 = tlefile.read("NOAA 20", str(fake_tlefile))

    assert tle_n20.line1 == "1 43013U 17073A   23045.54907786  .00000253  00000+0  14081-3 0  9995"
    assert tle_n20.line2 == "2 43013  98.7419 345.5839 0001610  80.3742 279.7616 14.19558274271576"


@pytest.mark.parametrize(("sat_name" ,"expected"),
                         [("NOAA 21",
                           "NOAA 21 (JPSS-2)"),
                          ("NOAA 2",
                           "NOAA 21 (JPSS-2)"),
                          ("NOAA 2",
                           "NOAA 20"),
                          ("NOAA",
                           "NOAA 21 (JPSS-2)"),
                          ("N",
                           "NOAA 21 (JPSS-2)")
                          ]
                         )
def test_read_tlefile_non_standard_platform_name_matching_start_of_name_in_tlefile(sat_name, expected,
                                                                                   caplog,
                                                                                   monkeypatch,
                                                                                   fake_platforms_txt_file,
                                                                                   fake_tlefile):
    """Test create a tle-object by reading tle data from file.

    Use non-standard naming matching only the beginning of what is in the TLE files.
    """
    from pyorbital import tlefile

    path_to_platforms_txt_file = fake_platforms_txt_file.parent
    monkeypatch.setenv("PYORBITAL_CONFIG_PATH", str(path_to_platforms_txt_file))

    with pytest.raises(KeyError) as exc_info:
        with caplog.at_level(logging.DEBUG):
            _ = tlefile.read(sat_name, str(fake_tlefile))

    assert f"Found a possible match: {expected}?" in caplog.text
    assert str(exc_info.value) == f'"Found no TLE entry for \'{sat_name}\'"'


@pytest.fixture
def fake_platforms_file(tmp_path):
    """Return file path to a fake platforms.txt file."""
    file_path = tmp_path / "platforms.txt"
    lines = ["# Some header lines - line 1\n",
             "# Some header lines - line 2\n",
             "NOAA-21 54234\n",
             "NOAA-20 43013\n",
             "UNKNOWN SATELLITE 99999\n"
             ]
    with open(file_path, "w") as fpt:
        fpt.writelines(lines)

    return file_path


@pytest.fixture(scope="session")
def fake_local_tles_dir(tmp_path_factory):
    """Make a list of fake tle files in a directory."""
    tle_dir = tmp_path_factory.mktemp("tle_files")
    file_path = tle_dir / "tle-202211180230.txt"
    file_path.touch()
    time.sleep(1)
    file_path = tle_dir / "tle-202211180430.txt"
    file_path.touch()
    time.sleep(1)
    file_path = tle_dir / "tle-202211180630.txt"
    file_path.touch()
    time.sleep(1)
    file_path = tle_dir / "tle-202211180830.txt"
    file_path.touch()

    return tle_dir


@pytest.fixture
def _mock_env_ppp_config_dir(monkeypatch):
    """Mock environment variable PPP_CONFIG_DIR."""
    monkeypatch.setenv("PPP_CONFIG_DIR", "/path/to/old/mpop/config/dir")


@pytest.fixture
def _mock_env_ppp_config_dir_missing(monkeypatch):
    """Mock that the environment variable PPP_CONFIG_DIR is missing."""
    monkeypatch.delenv("PPP_CONFIG_DIR", raising=False)


@pytest.fixture
def _mock_env_tles_missing(monkeypatch):
    """Mock that the environment variable TLES is missing."""
    monkeypatch.delenv("TLES", raising=False)


@pytest.fixture
def _mock_env_tles(monkeypatch, fake_local_tles_dir):
    """Mock environment variable TLES."""
    monkeypatch.setenv("TLES", os.path.join(fake_local_tles_dir, "*"))


@pytest.mark.usefixtures("_mock_env_ppp_config_dir_missing")
def test_get_config_path_no_env_defined(caplog):
    """Test getting the config path."""
    from pyorbital.tlefile import PKG_CONFIG_DIR, _get_config_path

    with caplog.at_level(logging.WARNING):
        res = _get_config_path()

    assert res == PKG_CONFIG_DIR
    assert caplog.text == ""


@pytest.mark.usefixtures("_mock_env_ppp_config_dir_missing")
def test_check_is_platform_supported_existing(caplog):
    """Test the function to check if an existing platform is supported on default."""
    from pyorbital.tlefile import PKG_CONFIG_DIR, check_is_platform_supported

    with caplog.at_level(logging.INFO):
        check_is_platform_supported("NOAA-21")

    logoutput_lines = caplog.text.split("\n")

    expected1 = "Satellite NOAA-21 is supported. NORAD number: 54234"
    expected2 = "Satellite names and NORAD numbers are defined in {path}".format(path=PKG_CONFIG_DIR)

    assert expected1 in logoutput_lines[0]
    assert expected2 in logoutput_lines[1]


@pytest.mark.usefixtures("_mock_env_ppp_config_dir_missing")
def test_check_is_platform_supported_unknown(caplog):
    """Test the function to check if an unknown  platform is supported on default."""
    from pyorbital.tlefile import PKG_CONFIG_DIR, check_is_platform_supported

    sat = "UNKNOWN"
    with caplog.at_level(logging.INFO):
        check_is_platform_supported(sat)

    logoutput_lines = caplog.text.split("\n")

    expected1 = "Satellite {satellite} is NOT supported.".format(satellite=sat)
    expected2 = ("Please add it to a local copy of the platforms.txt file and put in " +
                 "the directory pointed to by the environment variable PYORBITAL_CONFIG_PATH")
    expected3 = "Satellite names and NORAD numbers are defined in {path}".format(path=PKG_CONFIG_DIR)

    assert expected1 in logoutput_lines[0]
    assert expected2 in logoutput_lines[1]
    assert expected3 in logoutput_lines[2]


def test_get_config_path_ppp_config_set_but_not_pyorbital_future(caplog, monkeypatch):
    """Test getting the config path."""
    from pyorbital.tlefile import PKG_CONFIG_DIR, _get_config_path

    monkeypatch.setenv("SATPY_CONFIG_PATH", "/path/to/satpy/etc")
    monkeypatch.setenv("PPP_CONFIG_DIR", "/path/to/old/mpop/config/dir")

    with caplog.at_level(logging.WARNING):
        res = _get_config_path()

    log_output = ("The use of PPP_CONFIG_DIR is no longer supported! " +
                  "Please use PYORBITAL_CONFIG_PATH if you need a custom config path for pyorbital!")
    assert log_output in caplog.text
    assert res == PKG_CONFIG_DIR


def test_get_config_path_ppp_config_set_and_pyorbital(caplog, monkeypatch):
    """Test getting the config path."""
    from pyorbital.tlefile import _get_config_path

    pyorbital_config_dir = "/path/to/pyorbital/config/dir"
    monkeypatch.setenv("PYORBITAL_CONFIG_PATH", pyorbital_config_dir)
    monkeypatch.setenv("PPP_CONFIG_DIR", "/path/to/old/mpop/config/dir")

    with caplog.at_level(logging.WARNING):
        res = _get_config_path()

    assert res == pyorbital_config_dir
    assert caplog.text == ""


@pytest.mark.usefixtures("_mock_env_ppp_config_dir_missing")
def test_get_config_path_pyorbital_ppp_missing(caplog, monkeypatch):
    """Test getting the config path.

    The old mpop PPP_CONFIG_PATH is not set but the PYORBITAL one is.
    """
    from pyorbital.tlefile import _get_config_path

    pyorbital_config_dir = "/path/to/pyorbital/config/dir"
    monkeypatch.setenv("PYORBITAL_CONFIG_PATH", pyorbital_config_dir)

    with caplog.at_level(logging.DEBUG):
        res = _get_config_path()

    assert res == pyorbital_config_dir
    log_output = ("Path to the Pyorbital configuration (where e.g. " +
                  "platforms.txt is found): {path}".format(path=pyorbital_config_dir))
    assert log_output in caplog.text


def test_read_platform_numbers(fake_platforms_file):
    """Test reading the platform names and associated catalougue numbers."""
    from pyorbital.tlefile import read_platform_numbers

    res = read_platform_numbers(str(fake_platforms_file))
    assert res == {"NOAA-21": "54234", "NOAA-20": "43013", "UNKNOWN SATELLITE": "99999"}


@pytest.mark.usefixtures("_mock_env_tles_missing")
def test_get_local_tle_path_tle_env_missing():
    """Test getting the path to local TLE files - env TLES missing."""
    from pyorbital.tlefile import _get_local_tle_path_from_env

    res = _get_local_tle_path_from_env()
    assert res is None


@pytest.mark.usefixtures("_mock_env_tles")
def test_get_local_tle_path(fake_local_tles_dir):
    """Test getting the path to local TLE files."""
    from pyorbital.tlefile import _get_local_tle_path_from_env

    res = _get_local_tle_path_from_env()
    assert res == os.path.join(fake_local_tles_dir, "*")


def test_get_uris_and_open_func_using_tles_env(caplog, fake_local_tles_dir, monkeypatch):
    """Test getting the uris and associated open-function for reading tles.

    Test providing no tle file but using the TLES env to find local tle files.
    """
    from collections.abc import Sequence

    from pyorbital.tlefile import _get_uris_and_open_func

    monkeypatch.setenv("TLES", str(os.path.join(fake_local_tles_dir, "*")))

    with caplog.at_level(logging.DEBUG):
        uris, _ = _get_uris_and_open_func()

    assert isinstance(uris, Sequence)
    assert uris[0] == str(fake_local_tles_dir / "tle-202211180830.txt")
    log_message = "Reading TLE from {msg}".format(msg=str(fake_local_tles_dir))
    assert log_message in caplog.text


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
        assert tle.satnumber == "25544"
        assert tle.classification == "U"
        assert tle.id_launch_year == "98"
        assert tle.id_launch_number == "067"
        assert tle.id_launch_piece.strip() == "A"
        assert tle.epoch_year == "08"
        assert tle.epoch_day == 264.51782528
        epoch = (datetime.datetime(2008, 1, 1)
                 + datetime.timedelta(days=264.51782528 - 1))
        assert tle.epoch == epoch
        assert tle.mean_motion_derivative == -2.182e-05
        assert tle.mean_motion_sec_derivative == 0.0
        assert tle.bstar == -1.1606e-05
        assert tle.ephemeris_type == 0
        assert tle.element_number == 292

        # line 2
        assert tle.inclination == 51.6416
        assert tle.right_ascension == 247.4627
        assert tle.excentricity == 0.0006703
        assert tle.arg_perigee == 130.536
        assert tle.mean_anomaly == 325.0288
        assert tle.mean_motion == 15.72125391
        assert tle.orbit == 56353

    def test_from_line(self):
        """Test parsing from line elements."""
        from pyorbital.tlefile import Tle

        tle = Tle("ISS (ZARYA)", line1=LINE1, line2=LINE2)
        self.check_example(tle)

    def test_from_file(self):
        """Test reading and parsing from a file."""
        from pyorbital.tlefile import Tle

        filehandle, filename = mkstemp()
        try:
            os.write(filehandle, "\n".join([LINE0, LINE1, LINE2]).encode("utf-8"))
            os.close(filehandle)
            tle = Tle("ISS (ZARYA)", filename)
            self.check_example(tle)
        finally:
            os.remove(filename)

    def test_from_file_with_hyphenated_platform_name(self):
        """Test reading and parsing from a file with a slightly different name."""
        from pyorbital.tlefile import Tle

        filehandle, filename = mkstemp()
        try:
            os.write(filehandle, NOAA19_3LINES.encode("utf-8"))
            os.close(filehandle)
            tle = Tle("NOAA-19", filename)
            assert tle.satnumber == "33591"
        finally:
            os.remove(filename)

    def test_from_file_with_no_platform_name(self):
        """Test reading and parsing from a file with a slightly different name."""
        from pyorbital.tlefile import Tle

        filehandle, filename = mkstemp()
        try:
            os.write(filehandle, NOAA19_2LINES.encode("utf-8"))
            os.close(filehandle)
            tle = Tle("NOAA-19", filename)
            assert tle.satnumber == "33591"
        finally:
            os.remove(filename)

    def test_from_mmam_xml(self):
        """Test reading from an MMAM XML file."""
        from tempfile import TemporaryDirectory

        from pyorbital.tlefile import Tle

        save_dir = TemporaryDirectory()
        with save_dir:
            fname = os.path.join(save_dir.name, "20210420_Metop-B_ADMIN_MESSAGE_NO_127.xml")
            with open(fname, "w") as fid:
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

    @mock.patch("pyorbital.tlefile.requests")
    def test_fetch_plain_tle_not_configured(self, requests):
        """Test downloading and a TLE file from internet."""
        requests.get = mock.MagicMock()
        requests.get.return_value = _get_req_response(200)

        # Not configured
        self.dl.config["downloaders"] = {}
        res = self.dl.fetch_plain_tle()
        assert res == {}
        requests.get.assert_not_called()

    @mock.patch("pyorbital.tlefile.requests")
    def test_fetch_plain_tle_two_sources(self, requests):
        """Test downloading and a TLE file from internet."""
        requests.get = mock.MagicMock()
        requests.get.return_value = _get_req_response(200)

        # Two sources, one with multiple locations
        self.dl.config["downloaders"] = FETCH_PLAIN_TLE_CONFIG

        res = self.dl.fetch_plain_tle()
        assert "source_1" in res
        assert len(res["source_1"]) == 3
        assert res["source_1"][0].line1 == LINE1
        assert res["source_1"][0].line2 == LINE2
        assert "source_2" in res
        assert len(res["source_2"]) == 1
        assert mock.call("mocked_url_1", timeout=15) in requests.get.mock_calls
        assert len(requests.get.mock_calls) == 4

    @mock.patch("pyorbital.tlefile.requests")
    def test_fetch_plain_tle_server_is_a_teapot(self, requests):
        """Test downloading a TLE file from internet."""
        requests.get = mock.MagicMock()
        # No data returned because the server is a teapot
        requests.get.return_value = _get_req_response(418)

        # Two sources, one with multiple locations
        self.dl.config["downloaders"] = FETCH_PLAIN_TLE_CONFIG

        res = self.dl.fetch_plain_tle()
        # The sources are in the dict ...
        assert len(res) == 2
        # ... but there are no TLEs
        assert len(res["source_1"]) == 0
        assert len(res["source_2"]) == 0

        assert mock.call("mocked_url_1", timeout=15) in requests.get.mock_calls
        assert len(requests.get.mock_calls) == 4

    @mock.patch("pyorbital.tlefile.requests")
    def test_fetch_spacetrack_login_fails(self, requests):
        """Test downloading TLEs from space-track.org."""
        mock_post = mock.MagicMock()
        mock_session = mock.MagicMock()
        mock_session.post = mock_post
        requests.Session.return_value.__enter__.return_value = mock_session

        self.dl.config["platforms"] = {
            25544: "ISS"
        }
        self.dl.config["downloaders"] = FETCH_SPACETRACK_CONFIG

        # Login fails, because the server is a teapot
        mock_post.return_value.status_code = 418
        res = self.dl.fetch_spacetrack()
        # Empty list of TLEs is returned
        assert res == []
        # The login was anyway attempted
        mock_post.assert_called_with(
            "https://www.space-track.org/ajaxauth/login",
            data={"identity": "username", "password": "passw0rd"})

    @mock.patch("pyorbital.tlefile.requests")
    def test_fetch_spacetrack_get_fails(self, requests):
        """Test downloading TLEs from space-track.org."""
        mock_post = mock.MagicMock()
        mock_get = mock.MagicMock()
        mock_session = mock.MagicMock()
        mock_session.post = mock_post
        mock_session.get = mock_get
        requests.Session.return_value.__enter__.return_value = mock_session

        self.dl.config["platforms"] = {
            25544: "ISS"
        }
        self.dl.config["downloaders"] = FETCH_SPACETRACK_CONFIG

        # Login works, but something is wrong (teapot) when asking for data
        mock_post.return_value.status_code = 200
        mock_get.return_value.status_code = 418
        res = self.dl.fetch_spacetrack()
        assert res == []
        mock_get.assert_called_with("https://www.space-track.org/"
                                    "basicspacedata/query/class/tle_latest/"
                                    "ORDINAL/1/NORAD_CAT_ID/25544/format/tle")

    @mock.patch("pyorbital.tlefile.requests")
    def test_fetch_spacetrack_success(self, requests):
        """Test downloading TLEs from space-track.org."""
        mock_post = mock.MagicMock()
        mock_get = mock.MagicMock()
        mock_session = mock.MagicMock()
        mock_session.post = mock_post
        mock_session.get = mock_get
        requests.Session.return_value.__enter__.return_value = mock_session

        tle_text = "\n".join((LINE0, LINE1, LINE2))
        self.dl.config["platforms"] = {
            25544: "ISS"
        }
        self.dl.config["downloaders"] = FETCH_SPACETRACK_CONFIG

        # Login works and data is received
        mock_post.return_value.status_code = 200
        mock_get.return_value.status_code = 200
        mock_get.return_value.text = tle_text
        res = self.dl.fetch_spacetrack()
        assert len(res) == 1
        assert res[0].line1 == LINE1
        assert res[0].line2 == LINE2

    def test_read_tle_files(self):
        """Test reading TLE files from a file system."""
        from tempfile import TemporaryDirectory

        tle_text = "\n".join((LINE0, LINE1, LINE2))

        save_dir = TemporaryDirectory()
        with save_dir:
            fname = os.path.join(save_dir.name, "tle_20200129_1600.txt")
            with open(fname, "w") as fid:
                fid.write(tle_text)

            # Add a non-existent file, it shouldn't cause a crash
            nonexistent = os.path.join(save_dir.name, "not_here.txt")
            # Use a wildcard to collect files (passed to glob)
            starred_fname = os.path.join(save_dir.name, "tle*txt")
            self.dl.config["downloaders"] = {
                "read_tle_files": {
                    "paths": [fname, nonexistent, starred_fname]
                }
            }
            res = self.dl.read_tle_files()
        assert len(res) == 2
        assert res[0].line1 == LINE1
        assert res[0].line2 == LINE2

    def test_read_xml_admin_messages(self):
        """Test reading TLE files from a file system."""
        from tempfile import TemporaryDirectory

        save_dir = TemporaryDirectory()
        with save_dir:
            fname = os.path.join(save_dir.name, "20210420_Metop-B_ADMIN_MESSAGE_NO_127.xml")
            with open(fname, "w") as fid:
                fid.write(tle_xml)
            # Add a non-existent file, it shouldn't cause a crash
            nonexistent = os.path.join(save_dir.name, "not_here.txt")
            # Use a wildcard to collect files (passed to glob)
            starred_fname = os.path.join(save_dir.name, "*.xml")
            self.dl.config["downloaders"] = {
                "read_xml_admin_messages": {
                    "paths": [fname, nonexistent, starred_fname]
                }
            }
            res = self.dl.read_xml_admin_messages()

        # There are two sets of TLEs in the file.  And as the same file is
        # parsed twice, 4 TLE objects are returned
        assert len(res) == 4
        assert res[0].line1 == LINE1
        assert res[0].line2 == LINE2
        assert res[1].line1 == LINE1_2
        assert res[1].line2 == LINE2_2


def _get_req_response(code):
    req = mock.MagicMock()
    req.status_code = code
    req.text = "\n".join((LINE0, LINE1, LINE2))
    return req


class TestSQLiteTLE(unittest.TestCase):
    """Test saving TLE data to a SQLite database."""

    def setUp(self):
        """Create a database instance."""
        from tempfile import TemporaryDirectory

        from pyorbital.tlefile import SQLiteTLE, Tle

        self.temp_dir = TemporaryDirectory()
        self.db_fname = os.path.join(self.temp_dir.name, "tle.db")
        self.platforms = {25544: "ISS"}
        self.writer_config = {
            "output_dir": os.path.join(self.temp_dir.name, "tle_dir"),
            "filename_pattern": "tle_%Y%m%d_%H%M%S.%f.txt",
            "write_name": True,
            "write_always": False
        }
        self.db = SQLiteTLE(self.db_fname, self.platforms, self.writer_config)
        self.tle = Tle("ISS", line1=LINE1, line2=LINE2)

    def tearDown(self):
        """Clean temporary files."""
        with suppress(PermissionError, NotADirectoryError):
            self.temp_dir.cleanup()

    def test_init(self):
        """Test that the init did what it should have."""
        from pyorbital.tlefile import PLATFORM_NAMES_TABLE, table_exists

        columns = [col.strip() for col in
                   PLATFORM_NAMES_TABLE.strip("()").split(",")]
        num_columns = len(columns)

        assert os.path.exists(self.db_fname)
        assert table_exists(self.db.db, "platform_names")
        res = self.db.db.execute("select * from platform_names")
        names = [description[0] for description in res.description]
        assert len(names) == num_columns
        for col in columns:
            assert col.split(" ")[0] in names

    def test_update_db(self):
        """Test updating database with new data."""
        from pyorbital.tlefile import ISO_TIME_FORMAT, SATID_TABLE, _utcnow, table_exists

        # Get the column names
        columns = [col.strip() for col in
                   SATID_TABLE.replace("'{}' (", "").strip(")").split(",")]
        # Platform number
        satid = int(list(self.platforms.keys())[0])

        # Data from a platform that isn't configured
        self.db.platforms = {}
        self.db.update_db(self.tle, "foo")
        assert not table_exists(self.db.db, satid)
        assert not self.db.updated

        # Configured platform
        self.db.platforms = self.platforms
        self.db.update_db(self.tle, "foo")
        assert table_exists(self.db.db, satid)
        assert self.db.updated

        # Check that all the columns were added
        res = self.db.db.execute(f"select * from '{satid:d}'")  # noseq
        names = [description[0] for description in res.description]
        for col in columns:
            assert col.split(" ")[0] in names

        # Check the data
        data = res.fetchall()
        assert len(data) == 1
        # epoch
        assert data[0][0] == "2008-09-20T12:25:40.104192"
        # TLE
        assert data[0][1] == "\n".join((LINE1, LINE2))
        # Date when the data were added should be close to current time
        date_added = datetime.datetime.strptime(data[0][2], ISO_TIME_FORMAT)
        now = _utcnow()
        assert (now - date_added).total_seconds() < 1.0
        # Source of the data
        assert data[0][3] == "foo"

        # Try to add the same data again. Nothing should change even
        # if the source is different if the epoch is the same
        self.db.update_db(self.tle, "bar")
        res = self.db.db.execute(f"select * from '{satid:d}'")  # noseq
        data = res.fetchall()
        assert len(data) == 1
        date_added2 = datetime.datetime.strptime(data[0][2], ISO_TIME_FORMAT)
        assert date_added == date_added2
        # Source of the data
        assert data[0][3] == "foo"

    def test_write_tle_txt(self):
        """Test reading data from the database and writing it to a file."""
        import glob
        tle_dir = self.writer_config["output_dir"]

        # Put some data in the database
        self.db.update_db(self.tle, "foo")

        # Fake that the database hasn't been updated
        self.db.updated = False

        # Try to dump the data to disk
        self.db.write_tle_txt()

        # The output dir hasn't been created
        assert not os.path.exists(tle_dir)

        self.db.updated = True
        self.db.write_tle_txt()

        # The dir should be there
        assert os.path.exists(tle_dir)
        # There should be one file in the directory
        files = glob.glob(os.path.join(tle_dir, "tle_*txt"))
        assert len(files) == 1
        # The file should have been named with the date ('%' characters
        # not there anymore)
        assert "%" not in files[0]
        # The satellite name should be in the file
        with open(files[0], "r") as fid:
            data = fid.read().split("\n")
        assert len(data) == 3
        assert "ISS" in data[0]
        assert data[1] == LINE1
        assert data[2] == LINE2

        # Call the writing again, nothing should be written. In
        # real-life this assumes a re-run has been done without new
        # TLE data
        self.db.updated = False
        self.db.write_tle_txt()
        files = glob.glob(os.path.join(tle_dir, "tle_*txt"))
        assert len(files) == 1

        # Force writing with every call
        # Do not write the satellite name
        self.db.writer_config["write_always"] = True
        self.db.writer_config["write_name"] = False
        # Wait a bit to ensure different filename
        time.sleep(2)
        self.db.write_tle_txt()
        files = sorted(glob.glob(os.path.join(tle_dir, "tle_*txt")))
        assert len(files) == 2
        with open(files[1], "r") as fid:
            data = fid.read().split("\n")
        assert len(data) == 2
        assert data[0] == LINE1
        assert data[1] == LINE2

def test_tle_instance_printing():
    """Test the print the Tle instance."""
    from pyorbital.tlefile import Tle

    tle = Tle("ISS", line1=LINE1, line2=LINE2)

    expected = "{'arg_perigee': 130.536,\n 'bstar': -1.1606e-05,\n 'classification': 'U',\n 'element_number': 292,\n 'ephemeris_type': 0,\n 'epoch': np.datetime64('2008-09-20T12:25:40.104192'),\n 'epoch_day': 264.51782528,\n 'epoch_year': '08',\n 'excentricity': 0.0006703,\n 'id_launch_number': '067',\n 'id_launch_piece': 'A  ',\n 'id_launch_year': '98',\n 'inclination': 51.6416,\n 'mean_anomaly': 325.0288,\n 'mean_motion': 15.72125391,\n 'mean_motion_derivative': -2.182e-05,\n 'mean_motion_sec_derivative': 0.0,\n 'orbit': 56353,\n 'right_ascension': 247.4627,\n 'satnumber': '25544'}"  # noqa

    assert str(tle) == expected
