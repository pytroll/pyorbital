#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pytroll Community

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

import json
from datetime import timedelta
import numpy as np
from pyresample.spherical import SCoordinate
import pyresample as pr

# Norrköping coordinates:
NRK_LON = 16.1465
NRK_LAT = 58.5780
NRK_ALT = 0.03


def create_geojson_line(filename, arc):
    """From a pyresample arc vector store it to a Geojson file."""
    # geojson = {"type": "Feature"}
    geojson = {"type": "Feature", "geometry": {"type": "LineString",
                                               "coordinates": [
                                                   [arc.start.vertices_in_degrees[0][0],
                                                    arc.start.vertices_in_degrees[0][1]],
                                                   [arc.end.vertices_in_degrees[0][0],
                                                    arc.end.vertices_in_degrees[0][1]]
                                               ]}}

    with open(filename, 'w') as fp:
        json.dump(geojson, fp)


def get_arc_vector(timeobj, delta_t, sat, arc_len_min):
    """Get the arc defining the sub-satellite track in a certin time window.

    The time window is given by the start time 'timeobj' to start time
    'timeobj' plus time step 'delta_t'.
    """
    # Get positions at several points between +/- delta_t:
    len_one_arc_s = 60 * arc_len_min  # s
    num = int(delta_t.seconds/len_one_arc_s)
    # print(2*num)
    pos_vec = [sat.get_lonlatalt(timeobj + ind * 1.0/num * delta_t) for ind in range(-num, num + 1, 1)]

    # Calculate the Arc for each pixel. Later we sill see if arcs cross each other.
    # We could use only start and end. But then the SNO would be approximated some times to much.
    arcs = [pr.spherical.Arc(SCoordinate(lon=np.deg2rad(point1[0]),
                                         lat=np.deg2rad(point1[1])),
                             SCoordinate(lon=np.deg2rad(point2[0]),
                                         lat=np.deg2rad(point2[1])))
            for point1, point2 in zip(pos_vec[0:-1], pos_vec[1:])]

    return arcs


def get_sno_point(calipso, the_other_one, arc_calipso, arc_the_other_one, tobj, minthr):
    """Get the SNO point if there is any.

    If the two sub-satellite tracks of the overpasses intersects
    get the sub-satellite position and time where they cross,
    and determine if the time deviation is smaller than the require threshold:
    """
    import math
    intersect = arc_calipso.intersection(arc_the_other_one)
    point = (math.degrees(intersect.lon),
             math.degrees(intersect.lat))
    nextp = the_other_one.get_next_passes(tobj - timedelta(seconds=60*60),
                                          # SNO around tobj check for passes between +- one hour.
                                          # So that wanted pass for sure is next pass!
                                          2,  # Number of hours to find overpasses
                                          point[0],
                                          point[1],
                                          0)

    minthr_step = 20  # min less than half an orbit probably
    dtime = timedelta(seconds=60 * minthr_step * 2.0)

    if len(nextp) > 0:
        riset, fallt, maxt = nextp[0]
    else:
        print("No next passes found for, probably a bug!")
        tobj = tobj + dtime
        return None

    nextp = calipso.get_next_passes(tobj - timedelta(seconds=60*60),
                                    2,
                                    point[0],
                                    point[1],
                                    0)
    if len(nextp) > 0:
        riset, fallt, maxt_calipso = nextp[0]
    else:
        print("No next passes found for, probably a bug!")
        tobj = tobj + dtime
        return None

    # Get observer look from Norrkoping to the satellite when it is
    # in zenith over the SNO point:
    azi, elev = the_other_one.get_observer_look(maxt,
                                                NRK_LON,
                                                NRK_LAT,
                                                NRK_ALT)
    isNorrk = (elev > 0.0)

    tdelta = (maxt_calipso - maxt)
    tdmin = (tdelta.seconds + tdelta.days * 24*3600) / 60.

    if abs(tdmin) < minthr:
        match = {}
        match['satAdatetime'] = maxt
        match['satBdatetime'] = maxt_calipso
        match['sno_longitude'] = point[0]
        match['sno_latitude'] = point[1]
        match['minutes_diff'] = tdmin
        match['within_local_reception_area'] = isNorrk
        return match
