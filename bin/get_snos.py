#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pytroll

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

"""
Getting SNOs for two configurable satellite platforms.

SNO = Simultaneous Nadir Overpass: When two platforms sub-satellite track nadir
views cross each other in space and time. One can set a threshold in time
allowing the times to differ by up to a few minutes.

"""

import sys
import os

from datetime import datetime, timedelta

from pyorbital.sno_utils import get_arc_vector
from pyorbital.sno_utils import get_sno_point
from pyorbital.sno_utils import create_geojson_line


from pyorbital.orbital import Orbital
import numpy as np
import pandas as pd
from pyorbital.tlefile import Tle
import time
import logging


LOG = logging.getLogger('snos')
handler = logging.StreamHandler(sys.stderr)
handler.setLevel(0)
LOG.setLevel(0)
LOG.addHandler(handler)

tic = time.time()

R = 6370997.0

PATTERN = [0, 1, -1, 2, -2]
# TLE_DIR = '/data/orbital_elements/TLE/'
# TLE_DIR = '/home//temp'
TLE_DIR = "/home/a000680/data/tles"
TLE_SUBD_FORMAT = '%Y%m'
TLE_FILE_FORMAT = 'tle-%Y%m%d????.txt'
TLE_FILE_FORMAT2 = 'tle-%Y%m%d.txt'
TLE_FILE_FORMAT3 = 'TLE_noaa18.txt'


class NoTleFile(Exception):
    """Exception class to catch exceptions when a TLE file is not found."""

    pass


TLE_BUFFER_OTHER = {}
TLE_BUFFER_CALIPSO = {}
TLE_SATNAME = {'npp': 'SUOMI NPP',
               'snpp': 'SUOMI NPP',
               'aqua': 'AQUA',
               'metopb': 'METOP-B',
               'metopa': 'METOP-A',
               'Metop-C': 'METOP-C',
               'Metop-B': 'METOP-B',
               'Metop-A': 'METOP-A',
               'noaa19': 'NOAA 19',
               'noaa18': 'NOAA 18',
               'sentinel3a': 'SENTINEL-3A',
               'sentinel3b': 'SENTINEL-3B',
               'fengyun3d': 'FENGYUN 3D',
               'noaa15': 'NOAA 15',
               'noaa16': 'NOAA 16',
               'NOAA-18': 'NOAA 18',
               'NOAA-19': 'NOAA 19'
               }

TLE_IDs = {'npp': '37849',
           'snpp': '37849',
           'aqua': 'xxxxx',
           'calipso': '29108',
           'metopb': '38771',
           'Metop-B': '38771',
           'metopa': '29499',
           'noaa19': '33591',
           'noaa18': '28654',
           'sentinel3a': 'xxxxx',
           'sentinel3b': 'xxxxx',
           'fengyun3d': 'xxxxx',
           'noaa15': 'xxxxx',
           'noaa18': '28654',
           'noaa19': '33591',
           'NOAA-18': '28654',
           'NOAA-19': '33591',
           'noaa16': 'xxxxx'}

max_tle_days_diff = 3


def _get_tle_file(timestamp):
    # Find a not too old TLE file
    import glob
    for delta in PATTERN:
        dt_obj = timestamp - timedelta(days=delta)
        path = os.path.join(TLE_DIR, dt_obj.strftime(TLE_SUBD_FORMAT))
        if os.path.isdir(path):
            search = os.path.join(path, dt_obj.strftime(TLE_FILE_FORMAT))
            flst = glob.glob(search)
            if len(flst) != 0:
                fname = flst[0]
                LOG.info("Found TLE file: '%s'" % fname)
                return fname
            elif 1 == 2:
                search = os.path.join(path, dt_obj.strftime(TLE_FILE_FORMAT2))
                flst = glob.glob(search)
                if len(flst) != 0:
                    fname = flst[0]
                    LOG.info("Found TLE file: '%s'" % fname)
                    return fname

    raise NoTleFile("Found no TLE file close in time to " +
                    str(timestamp.strftime(TLE_FILE_FORMAT)))


def get_tle(platform, timestamp=None):
    """Get the tle from file, if not loaded already."""
    stamp = platform + timestamp.strftime('-%Y%m%d')
    try:
        tle = TLE_BUFFER_OTHER[stamp]
    except KeyError:
        tle = Tle(TLE_SATNAME[platform], _get_tle_file(timestamp))
        TLE_BUFFER_OTHER[stamp] = tle
    return tle


def get_tle_archive_calipso(timestamp, tle_filepath):
    """Get the TLEs from the archive of TLEs for the *calipso* platform."""
    return get_tle_archive(timestamp, tle_filepath, TLE_ID_CALIPSO, TLE_BUFFER_CALIPSO)


def get_tle_archive_other(timestamp, tle_filepath):
    """Get the TLEs from the archive of TLEs for the *other* platform."""
    return get_tle_archive(timestamp, tle_filepath, TLE_ID_OTHER, TLE_BUFFER_OTHER)


def populate_tle_buffer(filename, TLE_ID, MY_TLE_BUFFER):
    """Populate the TLE buffer."""
    with open(filename, 'r') as fh_:
        tle_data_as_list = fh_.readlines()
        for ind in range(0, len(tle_data_as_list), 2):
            if TLE_ID in tle_data_as_list[ind]:
                tle = Tle(TLE_ID, line1=tle_data_as_list[ind], line2=tle_data_as_list[ind+1])
                # dto = datetime.strptime(tle.epoch, '%Y-%m-%dT%H:%M:%S:%f')
                ts = (tle.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
                tobj = datetime.utcfromtimestamp(ts)
                MY_TLE_BUFFER[tobj] = tle


def get_tle_archive(timestamp, filename, TLE_ID, MY_TLE_BUFFER):
    """Get Two_Line elements from the archive.

    The TLE buffer MY_TLE_BUFFER is being updated.
    """
    # read tle data if not already in buffer
    if len(MY_TLE_BUFFER) == 0:
        populate_tle_buffer(filename, TLE_ID, MY_TLE_BUFFER)

    for tobj in MY_TLE_BUFFER:
        if tobj > timestamp:
            deltat = tobj - timestamp
        else:
            deltat = timestamp - tobj
        if np.abs((deltat).days) < 1:
            return MY_TLE_BUFFER[tobj]
    for delta_days in range(1, max_tle_days_diff + 1, 1):
        for tobj in MY_TLE_BUFFER:
            if tobj > timestamp:
                deltat = tobj - timestamp
            else:
                deltat = timestamp - tobj
            if np.abs((deltat).days) <= delta_days:
                print("Did not find TLE for {:s}, Using TLE from {:s}".format(tobj.strftime("%Y%m%d"),
                                                                              timestamp.strftime("%Y%m%d")))
                return MY_TLE_BUFFER[tobj]
    print("Did not find TLE for {:s} +/- 3 days")


def get_arguments():
    """Get the comman line arguments required to run the script."""
    import argparse

    parser = argparse.ArgumentParser(description='Calculate SNOS between two satellite platforms')
    parser.add_argument("-s", "--start-datetime",
                        required=True,
                        dest="start_datetime",
                        type=str,
                        default=None,
                        help="The datetime string corresponding to the start time of when SNOS should be calculated")
    parser.add_argument("-e", "--end-datetime",
                        required=True,
                        dest="end_datetime",
                        type=str,
                        default=None,
                        help="The datetime string corresponding to the end time of when SNOS should be calculated")
    parser.add_argument("-t", "--time-window",
                        required=True,
                        dest="time_window",
                        type=int,
                        default=None,
                        help=("The time window in number of minutes - the maximum time allowed between " +
                              "the two SNO observations"))
    parser.add_argument("--platform-name-A",
                        required=True,
                        dest="platform_name_one",
                        type=str,
                        default=None,
                        help="The name of the satellite platform (A)")
    parser.add_argument("--platform-name-B",
                        required=True,
                        dest="platform_name_two",
                        type=str,
                        default=None,
                        help="The name of the satellite platform (B)")
    parser.add_argument("--tolerance",
                        required=True,
                        dest="arc_len_min",
                        type=int,
                        default=None,
                        help=("The length in minutes of the delta arc used to find the point where "
                              "the two sub-satellite tracks cross. (2 is good, 5 is fast). "
                              "The lower this is the more accurate the SNO finding. The "
                              "higher this value is the less accurate but the faster the SNO finder is."))

    # approximation_arc_minutes 2-good 5-fast>"))

    # parser.add_argument("-c", "--configfile",
    #                     required=True,
    #                     dest="configfile",
    #                     type=str,
    #                     default=None,
    #                     help="The path to the configuration file")
    # parser.add_argument("-l", "--log-file", dest="log",
    #                     type=str,
    #                     default=None,
    #                     help="The file to log to (stdout per default).")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """Find SNOs for the two platforms within the time period given."""
    args = get_arguments()
    # platform_id_one = args.platform_name_one
    platform_id = args.platform_name_one
    # platform_id_two = args.platform_name_two
    calipso_id = args.platform_name_two

    minthr = int(args.time_window)
    time_start = datetime.strptime(args.start_datetime, "%Y%m%d%H%M")
    time_end = datetime.strptime(args.end_datetime, "%Y%m%d%H%M")
    arc_len_min = args.arc_len_min

    if platform_id not in TLE_SATNAME:
        LOG.error("Platform %s not supported!" % platform_id)
        sys.exit()

    # filename_calipso = "/home/a000680/data/tles/TLE_{:s}.txt".format(calipso_id)
    filename_calipso = "/home/a000680/data/tles/{:s}.tle".format(calipso_id)
    TLE_ID_CALIPSO = TLE_IDs.get(calipso_id)
    # filename_other = "/home/a000680/data/tles/TLE_{:s}.txt".format(platform_id)
    filename_other = "/home/a000680/data/tles/{:s}.tle".format(platform_id)
    TLE_ID_OTHER = TLE_IDs.get(platform_id)

    minthr_step = 20  # min less than half an orbit probably
    dtime = timedelta(seconds=60 * minthr_step * 2.0)
    # timestep_half = timedelta(seconds = 60 * minthr_step * 0.5)
    timestep_double = timedelta(seconds=60 * minthr_step * 2.0)
    # timestep = timedelta(seconds = 60 * minthr_step * 1.0)
    # make sure the two sat pass the SNO in the same step. We need and overlap of at least half minthr minutes.
    timestep_plus_30s = timedelta(seconds=60 * minthr_step * 1.0 + (minthr*0.5)*60 + 30)
    # delta_t = timedelta(seconds = 1200)
    # delta_t = timedelta(seconds = 60 * (minthr_step-1) * 2.0)

    calipso_obj = datetime(1970, 1, 1)
    other_obj = datetime(1970, 1, 1)

    tobj = time_start
    tle_calipso = None
    tle_the_other_one = None
    tobj_tmp = time_start
    isNorrk = False
    i = 0
    t_diff = timedelta(days=1)
    results = []
    while tobj < time_end:
        i = i + 1
        if i == 100:
            message = time.time() - tic, "seconds", dtime
            print(message)

        if not tle_calipso or calipso_obj - tobj > t_diff or tobj - calipso_obj > t_diff:
            tle_calipso = get_tle_archive_calipso(tobj, filename_calipso)
        if tle_calipso:
            ts = (tle_calipso.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
            calipso_obj = datetime.utcfromtimestamp(ts)
        if not tle_the_other_one or other_obj - tobj > t_diff or tobj - other_obj > t_diff:
            tle_the_other_one = get_tle_archive_other(tobj, filename_other)
        if tle_the_other_one:
            ts = (tle_the_other_one.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
            other_obj = datetime.utcfromtimestamp(ts)

        calipso = Orbital(calipso_id,
                          line1=tle_calipso.line1,
                          line2=tle_calipso.line2)
        the_other_one = Orbital(platform_id,
                                line1=tle_the_other_one.line1,
                                line2=tle_the_other_one.line2)

        got_intersection_acurate = False
        arc_calipso_vector = get_arc_vector(tobj, timestep_plus_30s, calipso, arc_len_min)
        arc_the_other_one_vector = get_arc_vector(tobj, timestep_plus_30s, the_other_one,  arc_len_min)
        # Approximate tracks with one arc each arc_len_min minutes.
        # For each pair of arcs check if they intersect.
        # There is atmost one intersection. Quit when we find it.
        for arc_calipso in arc_calipso_vector:
            for arc_the_other_one in arc_the_other_one_vector:
                if arc_calipso.intersects(arc_the_other_one):
                    got_intersection_acurate = True
                if got_intersection_acurate:
                    break
            if got_intersection_acurate:
                break

        if got_intersection_acurate:
            sno = get_sno_point(calipso, the_other_one,
                                arc_calipso, arc_the_other_one,
                                tobj, minthr)

            if sno:
                # For debugging:
                create_geojson_line('./calipso_arc_%d.geojson' % i, arc_calipso)
                results.append(sno)

                seconds_a = int(sno['satAdatetime'].strftime("%S")) + \
                    float(sno['satAdatetime'].strftime("%f"))/1000000.
                seconds_b = int(sno['satBdatetime'].strftime("%S")) + \
                    float(sno['satBdatetime'].strftime("%f"))/1000000.

                print("  " +
                      str(sno['satBdatetime'].strftime("%Y%m%d %H:%M")) +
                      "%5.1fs" % seconds_b + " "*5 +
                      str(sno['satAdatetime'].strftime("%Y%m%d %H:%M")) +
                      "%5.1fs" % seconds_a +
                      " "*6 + "(%7.2f, %7.2f)" % (sno['sno_latitude'], sno['sno_longitude']) +
                      "   " + "%4.1f min" % (sno['minutes_diff']) + "   " + str(sno['within_local_reception_area'])
                      )

        tobj = tobj + timestep_double
        if tobj - tobj_tmp > timedelta(days=1):
            tobj_tmp = tobj
            print(tobj_tmp.strftime("%Y-%m-%d"))
            LOG.debug(tobj_tmp.strftime("%Y-%m-%d"))

    dfresults = pd.DataFrame(results)

    print(dfresults)
    print(str(results[0]))

    geojson = {"type": "FeatureCollection", "features": []}

    # go through dataframe, append entries to geojson format
    for _, row in dfresults.iterrows():
        feature = {"type": "Feature", "geometry": {"type": "Point",
                                                   "coordinates": [row['sno_longitude'],
                                                                   row['sno_latitude']]},
                   "properties": {"datetime1": row['satAdatetime'].isoformat(),
                                  "datetime2": row['satBdatetime'].isoformat(),
                                  "tdiff_min": row['minutes_diff'],
                                  "within_area": row["within_local_reception_area"]}}
        geojson['features'].append(feature)

    import json
    with open('results.geojson', 'w') as fp:
        json.dump(geojson, fp)
