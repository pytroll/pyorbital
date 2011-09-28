#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011 SMHI

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""Test cases from the AIAA article.
"""

# TODO: right formal unit tests.

from pyorbital.orbital import Orbital, OrbitElements, _SGDP4
from pyorbital import tlefile
import numpy as np
from datetime import timedelta, datetime

class LineOrbital(Orbital):
    """Read TLE lines instead of file.
    """
    def __init__(self, satellite, line1, line2):
        satellite = satellite.upper()
        self.satellite_name = satellite
        self.tle = tlefile.read(satellite, line1=line1, line2=line2)
        self.orbit_elements = OrbitElements(self.tle)
        self._sgdp4 = _SGDP4(self.orbit_elements)


def get_results(satnumber, delay):
    """Get expected results from result file.
    """
    with open("aiaa_results") as f_2:
        line = f_2.readline()
        while(line):
            if line.endswith(" xx\n") and int(line[:-3]) == satnumber:
                line = f_2.readline()
                while(not line.startswith("%.8f"%delay)):
                    line = f_2.readline()
                sline = line.split()
                if delay == 0:
                    utc_time = None
                else:
                    utc_time = datetime.strptime(sline[-1], "%H:%M:%S.%f")
                    utc_time = utc_time.replace(year=int(sline[-4]),
                                                month=int(sline[-3]),
                                                day=int(sline[-2]))
                return (float(sline[1]),
                        float(sline[2]),
                        float(sline[3]),
                        float(sline[4]),
                        float(sline[5]),
                        float(sline[6]),
                        utc_time)
            line = f_2.readline()
                              
                

if __name__ == '__main__':
    with open("SGP4-VER.TLE") as f__:
        test_line = f__.readline()
        while(test_line):
            if test_line.startswith("#"):
                pass
            if test_line.startswith("1 "):
                line1 = test_line
            if test_line.startswith("2 "):
                line2 = test_line[:69]
                times = str.split(test_line[69:])
                times = np.arange(float(times[0]),
                                  float(times[1]) + 1,
                                  float(times[2]))
                try:
                    o = LineOrbital("unknown", line1, line2)
                except:
                    print "WARNING: skipping deep space computations"
                    test_line = f__.readline()
                    continue
                print "*" * 80
                print "satnumber:", o.tle.satnumber
                for delay in times:
                    try:
                        test_time = timedelta(minutes=delay) + o.tle.epoch
                        pos, vel = o.get_position(test_time, False)
                        res = get_results(int(o.tle.satnumber), float(delay))
                    except:
                        print "WARNING: TODO"
                        break

                    print "delay:" , delay, "minutes"
                    print "shown: pos_x(km) pos_y(km) pos_z(km) vel_x(km/s) vel_y(km/s) vel_z(km/s) time"
                    print "expected:"
                    print res[0], res[1], res[2], res[3], res[4], res[5], res[6]
                    print "got:"
                    print pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], test_time
                    print "difference (expected-got):"
                    if res[6] is not None:
                        print res[0] - pos[0], res[1] - pos[1], res[2] - pos[2], res[3] - vel[0], res[4] - vel[1], res[5] - vel[2], (res[6] - test_time).microseconds
                    else:
                        print res[0] - pos[0], res[1] - pos[1], res[2] - pos[2], res[3] - vel[0], res[4] - vel[1], res[5] - vel[2]
                    print "-" * 80
            test_line = f__.readline()
