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

"""Functions to support handling many archived TLEs."""


import numpy as np
import datetime as dt
from datetime import timezone
from pyorbital.tlefile import Tle

max_tle_days_diff = 3


def populate_tle_buffer(filename, tle_id, tle_buffer):
    """Populate the TLE buffer."""
    with open(filename, 'r') as fh_:
        tle_data_as_list = fh_.readlines()
        for ind in range(0, len(tle_data_as_list), 2):
            if tle_id in tle_data_as_list[ind]:
                tle = Tle(tle_id, line1=tle_data_as_list[ind], line2=tle_data_as_list[ind+1])
                # dto = datetime.strptime(tle.epoch, '%Y-%m-%dT%H:%M:%S:%f')
                ts = (tle.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
                # tobj = dt.datetime.utcfromtimestamp(ts)
                tobj = dt.datetime.fromtimestamp(ts, tz=timezone.utc)
                tle_buffer[tobj] = tle


def get_tle_archive(time_requested, filename, tle_id, tle_buffer):
    """Get Two-Line elements from the archive.

    The TLE buffer tle_buffer is being updated.
    """
    # Read tle data if not already in buffer
    if len(tle_buffer) == 0:
        populate_tle_buffer(filename, tle_id, tle_buffer)

    for tobj in tle_buffer:
        if tobj > time_requested:
            deltat = tobj - time_requested
        else:
            deltat = time_requested - tobj
        if np.abs((deltat).days) < 1:
            return tle_buffer[tobj]

    for delta_days in range(1, max_tle_days_diff + 1, 1):
        for tobj in tle_buffer:
            if tobj > time_requested:
                deltat = tobj - time_requested
            else:
                deltat = time_requested - tobj
            if np.abs((deltat).days) <= delta_days:
                print("Did not find TLE for {:s}, Using TLE from {:s}".format(tobj.strftime("%Y%m%d"),
                                                                              time_requested.strftime("%Y%m%d")))
                return tle_buffer[tobj]
    print("Did not find TLE for {:s} +/- 3 days")


def get_datetime_from_tle(tle_obj):
    """Get the datetime from a TLE object."""
    ts = (tle_obj.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    # return dt.datetime.utcfromtimestamp(ts)
    return dt.datetime.fromtimestamp(ts, tz=timezone.utc)
