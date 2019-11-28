#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 - 2019 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Nina Håkansson <nina.hakansson@smhi.se>
#   Erik Johansson <erik.johansson@smhi.se>

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

"""Getting SNOs for two configurable satellite platforms."""

import sys
from datetime import datetime, timedelta
import math
from pyorbital.orbital import Orbital
from pyorbital.sno_utils import get_config
from pyorbital.sno_utils import get_arc
from pyorbital.sno_utils import get_tle
import logging
import time

t = time.time()
LOG = logging.getLogger('snos')

handler = logging.StreamHandler(sys.stderr)
handler.setLevel(0)
LOG.setLevel(0)
LOG.addHandler(handler)


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
                        type=str,
                        default=None,
                        help=("The time window in number of minutes - the maximum time allowed between " +
                              "the two SNO observations"))
    parser.add_argument("-p", "--platform-name",
                        required=True,
                        dest="platform_name",
                        type=str,
                        default=None,
                        help="The name of the satellite platform")
    parser.add_argument("-a", "--area-configfile",
                        required=True,
                        dest="area_def",
                        type=str,
                        default=None,
                        help="The path to the area configuration file")
    parser.add_argument("-c", "--configfile",
                        required=True,
                        dest="configfile",
                        type=str,
                        default=None,
                        help="The path to the configuration file")
    parser.add_argument("-l", "--log-file", dest="log",
                        type=str,
                        default=None,
                        help="The file to log to (stdout per default).")

    args = parser.parse_args()
    return args


def get_sno_point(platform_id, cmp_platform, tobj, delta_t, arc_ref_pltfrm, arc_cmp_platform, station):
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

    nextp = ref_pltfrm.get_next_passes(tobj-delta_t,
                                       2,
                                       point[0],
                                       point[1],
                                       0)
    if len(nextp) > 0:
        riset, fallt, maxt_ref_pltfrm = nextp[0]
    else:
        LOG.warning("No next passes found! " + str(nextp))
        return None

    # Get observer look from the DR station to the satellite when it is
    # in at zenith over the SNO point:
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
    if abs(tdmin) < minthr:
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


if __name__ == "__main__":

    args = get_arguments()
    conf = get_config(args.configfile)

    station = {}
    station['lon'] = conf['station']['longitude']
    station['lat'] = conf['station']['latitude']
    station['alt'] = conf['station']['altitude']

    ref_platform_name = conf['reference_platform']['name']

    tle_dirs = conf['tle-dirs']
    tle_file_format = conf['tle-file-format']
    platform_id = args.platform_name
    minthr = int(args.time_window)
    time_start = datetime.strptime(args.start_datetime, "%Y%m%d")
    time_end = datetime.strptime(args.end_datetime, "%Y%m%d")

    dtime = timedelta(seconds=60 * minthr * 2.0)
    timestep = timedelta(seconds=60 * minthr * 1.0)
    delta_t = dtime

    tobj = time_start
    tle_ref_pltfrm = None
    tle_cmp_platform = None
    tobj_tmp = time_start
    is_within_antenna_horizon = False
    while tobj < time_end:
        if not tle_ref_pltfrm or abs(tle_ref_pltfrm.epoch.astype(datetime) - tobj) > timedelta(days=1):
            tle_ref_pltfrm = get_tle(tle_dirs, tle_file_format, ref_platform_name, tobj)
        if (not tle_cmp_platform or
                abs(tle_cmp_platform.epoch.astype(datetime) - tobj) > timedelta(days=2)):
            tle_cmp_platform = get_tle(tle_dirs, tle_file_format, platform_id, tobj)

        ref_pltfrm = Orbital(ref_platform_name,
                             line1=tle_ref_pltfrm.line1,
                             line2=tle_ref_pltfrm.line2)
        cmp_platform = Orbital(platform_id,
                               line1=tle_cmp_platform.line1,
                               line2=tle_cmp_platform.line2)

        arc_ref_pltfrm = get_arc(tobj, timestep, ref_pltfrm)
        arc_cmp_platform = get_arc(tobj, timestep, cmp_platform)

        if arc_ref_pltfrm and arc_cmp_platform:
            if arc_ref_pltfrm.intersects(arc_cmp_platform):
                # If the two sub-satellite tracks of the overpasses intersects
                # get the sub-satellite position and time where they cross,
                # and determine if the time deviation is smaller than the require threshold:
                sno = get_sno_point(platform_id, cmp_platform, tobj, delta_t, arc_ref_pltfrm,
                                    arc_cmp_platform, station)

                if sno:
                    print("  " +
                          str(sno['maxt_ref_pltfrm'].strftime("%Y-%m-%d %H:%M ")) +
                          "%5.1f" % sno['ref_pltfrmsec'] + " "*5 +
                          str(sno['maxt'].strftime("%Y-%m-%d %H:%M ")) +
                          "%5.1f" % sno['sec'] +
                          " imager-orbit %d " % (sno['orbit_nr'] + 1) +
                          " "*6 + "%7.2f  %7.2f" % (sno['point'][1], sno['point'][0]) +
                          "   " + "%4.1f" % abs(sno['tdmin']) + "   " + str(sno['is_within_antenna_horizon'])
                          )

        else:
            LOG.error("Failed getting the track-segments")

        tobj = tobj + dtime
        if tobj - tobj_tmp > timedelta(days=1):
            tobj_tmp = tobj
            LOG.debug(tobj_tmp.strftime("%Y-%m-%d"))
