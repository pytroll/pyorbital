#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pyorbital developers


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

"""Testing the TLE archive functions."""


import numpy as np
import datetime as dt
from datetime import timezone
import pyorbital
# from pyorbital.tle_archive import populate_tle_buffer
# from pyorbital.tle_archive import get_tle_archive
from pyorbital.tle_archive import TwoLineElementsFinder


def test_populate_tle_buffer(fake_tle_file1_calipso):
    """Test populate the TLE buffer."""
    tle_filename = str(fake_tle_file1_calipso)
    tleid_calipso = '29108'
    tle_finder = TwoLineElementsFinder(tleid_calipso, tle_filename)
    tle_finder.populate_tle_buffer()
    tlebuff = tle_finder.tle_buffer

    with open(tle_filename, 'r') as fpt:
        tlelines = fpt.readlines()

    expected_dtimes = [dt.datetime(2013, 12, 31, 13, 38, 47, 751936, tzinfo=timezone.utc),
                       dt.datetime(2014, 1, 1, 17, 39, 47, 34720, tzinfo=timezone.utc),
                       dt.datetime(2014, 1, 2, 20, 1, 53, 254560, tzinfo=timezone.utc),
                       dt.datetime(2014, 1, 4, 0, 2, 52, 226304, tzinfo=timezone.utc)]
    for idx, key in enumerate(tlebuff.keys()):
        assert key == expected_dtimes[idx]

    for idx, key in enumerate(tlebuff.keys()):
        tleobj = tlebuff[key]
        assert isinstance(tleobj, pyorbital.tlefile.Tle)
        assert tlelines[idx*2].strip() == tleobj.line1
        assert tlelines[idx*2+1].strip() == tleobj.line2


def test_get_tle_archive(fake_tle_file1_calipso):
    """Test getting all the TLEs from file with many TLEs."""
    tle_filename = str(fake_tle_file1_calipso)
    tleid_calipso = '29108'
    dtobj = dt.datetime(2014, 1, 2, tzinfo=timezone.utc)

    tle_finder = TwoLineElementsFinder(tleid_calipso, tle_filename)
    tleobj = tle_finder.get_tle_archive(dtobj)

    ts = (tleobj.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    dtime_valid = dt.datetime.fromtimestamp(ts, tz=timezone.utc)
    expected = dt.datetime(2014, 1, 1, 17, 39, 47, 34720, tzinfo=timezone.utc)

    assert abs(expected - dtime_valid).total_seconds() < 0.001

    dtobj = dt.datetime(2014, 1, 1, 12, tzinfo=timezone.utc)
    tle_finder = TwoLineElementsFinder(tleid_calipso, tle_filename)
    tleobj = tle_finder.get_tle_archive(dtobj)

    ts = (tleobj.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    dtime_valid = dt.datetime.fromtimestamp(ts, tz=timezone.utc)
    expected = dt.datetime(2013, 12, 31, 13, 38, 47, 751936, tzinfo=timezone.utc)

    assert abs(expected - dtime_valid).total_seconds() < 0.001

    dtobj = dt.datetime(2014, 1, 1, 16, tzinfo=timezone.utc)
    tle_finder = TwoLineElementsFinder(tleid_calipso, tle_filename)
    tleobj = tle_finder.get_tle_archive(dtobj)

    ts = (tleobj.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    dtime_valid = dt.datetime.fromtimestamp(ts, tz=timezone.utc)
    expected = dt.datetime(2014, 1, 1, 17, 39, 47, 34720, tzinfo=timezone.utc)

    assert abs(expected - dtime_valid).total_seconds() < 0.001
