#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013 Martin Raspaud

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

"""Simple usage for geoloc.
"""

import numpy as np
from datetime import datetime
from pyorbital.geoloc import ScanGeometry, compute_pixels, get_lonlatalt
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

# Couple of example Two Line Elements
tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"

# Choosing a specific time, this should be relatively close to the issue date of the TLE
t = datetime(2012, 12, 12, 4, 16, 1, 575000)
# this is the number of full scan rotations
scans_nb = 10
# we take only every 40th point for plotting clarity
scan_points = np.arange(24, 2048, 40)
# This the maximum scan angle away from nadir for the given TLE that still sees earth.
scan_angle = 55.37
# period of one full rotation 1/6 s
scan_p = 0.16666667
# integration time of instrument
int_t = 0.000025

# build the avhrr instrument (scan angles)
# creates list of radian angles centered around nadir based on the scan points that see earth
avhrr = np.vstack(((scan_points / 1023.5-1) * np.deg2rad(-scan_angle),
                   np.zeros((len(scan_points),))))
avhrr = np.tile(
        avhrr[:, np.newaxis, :], [1, scans_nb, 1])

# building the corresponding times array
times = np.tile(scan_points * int_t, [scans_nb, 1])
offset = np.arange(scans_nb) * scan_p
times += np.expand_dims(offset, 1)

# build the scan geometry object
sgeom = ScanGeometry(avhrr, times)

# roll, pitch, yaw in radians. This is a static offset.
rpy = (0, 0, 0)

# print the longitude and latitude for the pixel positions
s_times = sgeom.times(t)
pixels_pos = compute_pixels((tle1, tle2), sgeom, s_times, rpy)
pos_time = get_lonlatalt(pixels_pos, s_times)

print(pos_time)

# Plot the result
m = Basemap(projection='stere', llcrnrlat=24, urcrnrlat=70, llcrnrlon=-25, urcrnrlon=120,
            lat_ts=58, lat_0=58, lon_0=14, resolution='l')

# convert and plot the predicted pixels in red
x, y = m(pos_time[0], pos_time[1])
p1 = m.plot(x, y, marker='+', color='red', markerfacecolor='red', markeredgecolor='red', markersize=1, markevery=1,
            zorder=4, linewidth=0.0)
m.fillcontinents(color='0.85', lake_color=None, zorder=3)
m.drawparallels(np.arange(-90., 90., 5.), labels=[1, 0, 1, 0], fontsize=10, dashes=[1, 0],
                color=[0.8, 0.8, 0.8], zorder=1)
m.drawmeridians(np.arange(-180., 180., 5.), labels=[0, 1, 0, 1], fontsize=10, dashes=[1, 0],
                color=[0.8, 0.8, 0.8], zorder=2)

plt.show()
