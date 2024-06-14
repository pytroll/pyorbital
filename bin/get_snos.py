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

from pyorbital.sno_utils import SNOfinder
from pyorbital.logger import setup_logging_from_config
import datetime as dt
from datetime import timezone
import logging


logger = logging.getLogger('snos')

# LOG = logging.getLogger('snos')
# handler = logging.StreamHandler(sys.stderr)
# handler.setLevel(0)
# LOG.setLevel(0)
# LOG.addHandler(handler)


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
    parser.add_argument("-c", "--configfile",
                        required=True,
                        dest="configfile",
                        type=str,
                        default=None,
                        help="The path to the configuration file")
    parser.add_argument("-l", "--log-config", dest="log_config",
                        type=str,
                        default=None,
                        help="Log config file to use instead of the standard logging.")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """Find SNOs for the two platforms within the time period given."""
    args = get_arguments()
    if args.log_config:
        setup_logging_from_config(args)

    logger.info("Starting up.")

    platform_id_one = args.platform_name_one
    platform_id_two = args.platform_name_two

    minutes_thr = int(args.time_window)
    starttime = dt.datetime.strptime(args.start_datetime, "%Y%m%d%H%M")
    starttime = starttime.replace(tzinfo=timezone.utc)

    endtime = dt.datetime.strptime(args.end_datetime, "%Y%m%d%H%M")
    endtime = endtime.replace(tzinfo=timezone.utc)
    arclength_minutes = args.arc_len_min

    sno_finder = SNOfinder(platform_id_one, platform_id_two, (starttime, endtime),
                           minutes_thr, arclength_minutes)
    sno_finder.set_configuration(args.configfile)
    results = sno_finder.get_snos_within_time_window()

    logger.info("Finished getting SNOs")
    print(results)

    sno_finder.dataframe2geojson(results)
    sno_finder.write_geojson('./results.geojson')
