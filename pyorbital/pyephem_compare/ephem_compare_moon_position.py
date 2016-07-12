#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015, 2016 Adam.Dybbroe

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


"""Compare my derivation of the altitude and azimuth of the moon with that
derived from ephem

"""

import ephem
import datetime
import numpy as np
from pyorbital import planets

if __name__ == "__main__":

    # Norrköping: Lat N 58° 35′ 15″ Lon E 16° 11′ 15″
    llon = 16 + 11.0 / 60. + 15. / 3600.
    llat = 58. + 35. / 60. + 15. / 3600.

    currtime = datetime.datetime(2012, 1, 7, 12, 0)
    endtime = datetime.datetime(2012, 1, 8, 12, 0)
    delta_t = datetime.timedelta(seconds=600)
    heights = []
    dtimes = []

    norrk = ephem.Observer()
    norrk.lat, norrk.lon = '58.5875', '16.1875'
    norrk.pressure = 0
    #norrk.horizon = '-0:34'

    dtimes = []
    ephem_heights = []
    heights = []
    ephem_azimuths = []
    azimuths = []
    while currtime < endtime:
        norrk.date = currtime.strftime('%Y/%m/%d %H:%M')
        m = ephem.Moon(norrk)
        ephem_heights.append(np.rad2deg(m.alt))
        ephem_azimuths.append(np.rad2deg(m.az))

        moon = planets.Moon(currtime)
        rasc, decl, alt, azi = moon.topocentric_position(llon, llat)
        # azi, alt, ha_ = planets.azimuthal_coord(currtime,
        #                                        moon.lonsun,
        #                                        llon, llat, rasc, decl)
        heights.append(alt)
        azimuths.append(azi)
        dtimes.append(currtime)
        currtime = currtime + delta_t

    heights = np.array(heights)
    ephem_heights = np.array(ephem_heights)
    azimuths = np.array(azimuths)
    ephem_azimuths = np.array(ephem_azimuths)

    # Plot the moon positions (altitude above horizon):
    import matplotlib.pyplot as plt

    fig = plt.figure()

    ax = fig.add_subplot(111)
    ax.scatter(ephem_heights, ephem_heights - heights)

    ax.set_xlabel('Ephem moon height (deg)', fontsize=20)
    ax.set_ylabel('Difference - Ephem-pyorbital (deg)', fontsize=20)
    ax.set_title('Comparing Moon height derivation')
    ax.grid(True)

    plt.savefig('./moonheight_compare.png')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ephem_azimuths, ephem_azimuths - azimuths)

    ax.set_xlabel('Ephem moon azimuth (deg)', fontsize=20)
    ax.set_ylabel('Difference - Ephem-pyorbital (deg)', fontsize=20)
    ax.set_title('Comparing Moon azimuth derivation')
    ax.grid(True)

    plt.savefig('./moonazimuth_compare.png')
