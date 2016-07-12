#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>

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

"""Compare pyorbital and ephem to data found at mooncalc.org 
"""

import ephem
from datetime import datetime
from pyorbital.moon_phase import moon_phase
from pyorbital import planets

# lon = 16.183
# lat = 58.6

# From mooncalc.org:
# Moon rising:	11:40
# Moon peak:	17:07
# Moon set:	22:47
# Moon distance:	367697km
# Moon elevation angle:	-21.47°
# Moon horizontal angle:	288.28°
# Phase = 31.6%

lon = 10.39307
lat = 57.07657

# Moon elevation angle:	-23.09°
# Moon horizontal angle:	289.98°
# Phase = 31.9%

#currtime = datetime(2015, 12, 16, 23, 16)
currtime = datetime(2015, 12, 16, 23, 50)

# Pyorbital:
phase = moon_phase(currtime)
moon = planets.Moon(currtime)
rasc, decl, alt, azi = moon.topocentric_position(lon, lat)
print("pyorbital: phase=%f alt=%f azi=%f" % (phase * 100, alt, azi))

# Ephem:
norrk = ephem.Observer()
# norrk.lat, norrk.lon = '58.6', '16.183'
norrk.lat, norrk.lon = '57.07657', '10.39307'
norrk.pressure = 0
norrk.date = currtime.strftime('%Y/%m/%d %H:%M')
m = ephem.Moon(norrk)
print("ephem: phase=%f alt=%s azi=%s" %
      (m.phase, str(m.alt), str(m.az)))

# pyorbital: phase=31.638283 alt=-21.904673 azi=288.377773
# ephem: phase=31.796560 alt=-21:32:21.0 azi=288:19:25.5

# pyorbital: phase=31.884964 alt=-23.521642 azi=290.057501
# ephem: phase=32.042305 alt=-23:09:50.4 azi=290:01:13.8
