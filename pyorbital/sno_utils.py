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
import math
import logging
import yaml
from yaml import SafeLoader

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
    """Get the configuration from file.

    :configfile: The file path of the yaml configuration file.

    :return: A configuration dictionary.

    """
    config = {}
    with open(configfile, 'r') as fp_:
        config = _recursive_dict_update(config, yaml.load(fp_, Loader=SafeLoader))

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
    """Get the tle from file, if not loaded already.

    :tle_dirs: A list of paths where tle-files are located
    :tle_file_format: The tle-file format pattern
    :platorm: Satellite platform name

    :return: A pyorbital Tle object

    """
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

    :timeobj: Start time (datetime object)
    :delta_t: Time step to add to the start time
    :sat: Orbital object for the satellite platform

    :return: A Pyresample spherical geometry arc object

    """
    # Get the start and end positions:
    pos_start = sat.get_lonlatalt(timeobj - delta_t)
    pos_end = sat.get_lonlatalt(timeobj + delta_t)

    coord_start = Coordinate(lon=pos_start[0], lat=pos_start[1])
    coord_end = Coordinate(lon=pos_end[0], lat=pos_end[1])

    return Arc(coord_start, coord_end)


def get_sno_point(platform_id, cmp_platform, ref_platform, tobj, delta_t,
                  arc_ref_pltfrm, arc_cmp_platform, station, minute_thr):
    """Get the SNO point if there is any.

    If the two sub-satellite tracks of the overpasses intersects
    get the sub-satellite position and time where they cross,
    and determine if the time deviation is smaller than the require threshold:
    """
    intersect = arc_ref_pltfrm.intersection(arc_cmp_platform)
    point = (math.degrees(intersect.lon),
             math.degrees(intersect.lat))

    nextp = cmp_platform.get_next_passes(tobj-delta_t,
                                         2,
                                         point[0],
                                         point[1],
                                         0)
    if len(nextp) > 0:
        riset, fallt, maxt = nextp[0]
    else:
        LOG.warning("No next passes found for " +
                    platform_id + "! " + str(nextp))
        return None

    nextp = ref_platform.get_next_passes(tobj-delta_t,
                                         2,
                                         point[0],
                                         point[1],
                                         0)
    if len(nextp) > 0:
        riset, fallt, maxt_ref_pltfrm = nextp[0]
    else:
        LOG.warning("No next passes found! " + str(nextp))
        return None

    # Get observer look from the specified DR station to the satellite when it
    # is at zenith over the SNO point:
    azi, elev = cmp_platform.get_observer_look(maxt,
                                               station['lon'],
                                               station['lat'],
                                               station['alt'])

    is_within_antenna_horizon = (elev > 0.0)

    ref_pltfrmsec = (int(maxt_ref_pltfrm.strftime("%S")) +
                     int(maxt_ref_pltfrm.strftime("%f"))/1000000.)
    sec = (int(maxt.strftime("%S")) +
           int(maxt.strftime("%f"))/1000000.)

    tdelta = (maxt_ref_pltfrm - maxt)
    tdmin = (tdelta.seconds + tdelta.days * 24*3600) / 60.
    if abs(tdmin) < minute_thr:
        orbit_nr = cmp_platform.get_orbit_number(maxt)

        result = {}
        result['maxt'] = maxt
        result['maxt_ref_pltfrm'] = maxt_ref_pltfrm
        result['ref_pltfrmsec'] = ref_pltfrmsec
        result['sec'] = sec
        result['is_within_antenna_horizon'] = is_within_antenna_horizon
        result['tdmin'] = tdmin
        result['orbit_nr'] = orbit_nr
        result['point'] = point

        return result
    else:
        return None
