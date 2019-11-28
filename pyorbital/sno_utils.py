#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Pytroll

# Author(s):

#   Adam.Dybbroe <dam.dybbroe@smhi.se>

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

"""Utility functions to find simultaneous nadir overpasses (SNOs)."""

import os
from datetime import timedelta
from pyorbital.tlefile import Tle
from pyresample.spherical_geometry import Coordinate, Arc

import logging
import yaml
try:
    from yaml import UnsafeLoader
except ImportError:
    from yaml import Loader as UnsafeLoader
try:
    # python 3.3+
    from collections.abc import Mapping
except ImportError:
    # deprecated (above can't be done in 2.7)
    from collections import Mapping

LOG = logging.getLogger(__name__)

TLE_BUFFER = {}
TLE_SATNAME = {'Suomi-NPP': 'SUOMI NPP',
               'NOAA-20': 'NOAA 20',
               'EOS-Aqua': 'AQUA',
               'Metop-B': 'METOP-B',
               'NOAA-19': 'NOAA 19',
               'NOAA-18': 'NOAA 18'}


class NoTleFile(Exception):
    """Exception to catch missing TLE file."""

    pass


def _recursive_dict_update(d, u):
    """Recursive dictionary update.

    Copied from:

        http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    """
    for k, v in u.items():
        if isinstance(v, Mapping):
            r = _recursive_dict_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


def get_config(configfile):
    """Get the configuration from file."""
    config = {}
    with open(configfile, 'r') as fp_:
        config = _recursive_dict_update(config, yaml.load(fp_, Loader=UnsafeLoader))

    return config


def _get_tle_file(tledirs, tle_file_format, timestamp):
    # Find a not too old TLE file
    for path in tledirs:
        if os.path.isdir(path):
            for idx in range(5):
                dt_obj = timestamp - timedelta(days=idx)
                fname = os.path.join(path, dt_obj.strftime(tle_file_format))
                if os.path.isfile(fname):
                    LOG.info("Found TLE file: '%s'" % fname)
                    return fname
    raise NoTleFile("Found no TLE file close in time to " +
                    str(timestamp.strftime(tle_file_format)))


def get_tle(tle_dirs, tle_file_format, platform, timestamp=None):
    """Get the tle from file, if not loaded already."""
    stamp = platform + timestamp.strftime('-%Y%m%d')
    try:
        tle = TLE_BUFFER[stamp]
    except KeyError:
        tle = Tle(TLE_SATNAME.get(platform, platform),
                  _get_tle_file(tle_dirs, tle_file_format, timestamp))
        TLE_BUFFER[stamp] = tle
    return tle


def get_arc(timeobj, delta_t, sat):
    """Get the arc on the geoid of the sub-satellit track over a small time interval.

    It get's the arc defining the sub-satellite track from start time 'timeobj'
    to start time 'timeobj' plus a time step 'delta_t'.
    """
    # Get the start and end positions:
    pos_start = sat.get_lonlatalt(timeobj - delta_t)
    pos_end = sat.get_lonlatalt(timeobj + delta_t)

    coord_start = Coordinate(lon=pos_start[0], lat=pos_start[1])
    coord_end = Coordinate(lon=pos_end[0], lat=pos_end[1])

    return Arc(coord_start, coord_end)
