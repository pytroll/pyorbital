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
from __future__ import with_statement

import os

from pyorbital.orbital import Orbital, OrbitElements, _SGDP4
from pyorbital import tlefile, astronomy
import numpy as np
from datetime import timedelta, datetime
import unittest

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
    path = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(path, "aiaa_results")) as f_2:
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

class AIAAIntegrationTest(unittest.TestCase):
    """Test against the AIAA test cases.
    """
    
    def test_aiaa(self):
        """Do the tests against AIAA test cases.
        """
        path = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(path, "SGP4-VER.TLE")) as f__:
            test_line = f__.readline()
            while(test_line):
                if test_line.startswith("#"):
                    test_name = test_line
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
                    except Exception, e:
                        # WARNING: skipping deep space computations
                        from warnings import warn
                        warn(test_name + ' ' + str(e))
                        
                        test_line = f__.readline()
                        continue
                    for delay in times:
                        try:
                            test_time = timedelta(minutes=delay) + o.tle.epoch
                            pos, vel = o.get_position(test_time, False)
                            res = get_results(int(o.tle.satnumber), float(delay))
                        except (NotImplementedError, ValueError), e:
                            # WARNING: TODO
                            from warnings import warn
                            warn(test_name + ' ' + str(e))
                            break

                        delta_pos = 5e-6 # km =  5 mm
                        delta_vel = 5e-9 # km/s = 5 um/s
                        delta_time = 1e-3 # 1 milisecond
                        self.assertTrue(abs(res[0] - pos[0]) < delta_pos)
                        self.assertTrue(abs(res[1] - pos[1]) < delta_pos)
                        self.assertTrue(abs(res[2] - pos[2]) < delta_pos)
                        self.assertTrue(abs(res[3] - vel[0]) < delta_vel)
                        self.assertTrue(abs(res[4] - vel[1]) < delta_vel)
                        self.assertTrue(abs(res[5] - vel[2]) < delta_vel)
                        if res[6] is not None:
                            dt = astronomy._days(res[6] - test_time) * 24 * 60
                            self.assertTrue(abs(dt) < delta_time)
                        
                test_line = f__.readline()
        

if __name__ == '__main__':
    unittest.main()
