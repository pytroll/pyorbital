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


import datetime as dt
import pyorbital
from pyorbital.tle_archive import populate_tle_buffer


def test_populate_tle_buffer(fake_tle_file1_calipso):
    """Test populate the TLE buffer."""
    tle_filename = str(fake_tle_file1_calipso)
    tlebuff = {}
    tleid_calipso = '29108'
    populate_tle_buffer(tle_filename, tleid_calipso, tlebuff)

    with open(tle_filename, 'r') as fpt:
        tlelines = fpt.readlines()

    expected_dtimes = [dt.datetime(2013, 12, 31, 13, 38, 47, 751936),
                       dt.datetime(2014, 1, 1, 17, 39, 47, 34720),
                       dt.datetime(2014, 1, 2, 20, 1, 53, 254560),
                       dt.datetime(2014, 1, 4, 0, 2, 52, 226304)]
    for idx, key in enumerate(tlebuff.keys()):
        assert key == expected_dtimes[idx]

    for idx, key in enumerate(tlebuff.keys()):
        tleobj = tlebuff[key]
        assert isinstance(tleobj, pyorbital.tlefile.Tle)
        assert tlelines[idx*2].strip() == tleobj.line1
        assert tlelines[idx*2+1].strip() == tleobj.line2
