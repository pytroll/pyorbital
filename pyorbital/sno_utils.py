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
import datetime as dt
from datetime import timedelta, timezone
from geopy import distance
# from geojson import dump
import numpy as np
import pandas as pd

from pyresample.spherical import SCoordinate
import pyresample as pr
from pyorbital.config import get_config
from pyorbital.orbital import Orbital
from pyorbital.tlefile import SATELLITES
from pyorbital.tle_archive import get_tle_archive
from pyorbital.tle_archive import get_datetime_from_tle

from trollsift.parser import Parser
from pathlib import Path
import time
import logging


OSCAR_NAMES = {'npp': 'Suomi-NPP',
               'snpp': 'Suomi-NPP',
               'aqua': 'EOS-Aqua',
               'metopb': 'Metop-B',
               'metopa': 'Metop-A',
               'noaa19': 'NOAA-19',
               'noaa18': 'NOAA-18',
               'sentinel3a': 'Sentinel-3A',
               'sentinel3b': 'Sentinel-3B',
               'fengyun3d': 'FY-3D',
               'noaa15': 'NOAA-15',
               'noaa16': 'NOAA-16',
               'calipso': 'CALIPSO'
               }


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

ZERO_SECONDS = timedelta(seconds=0)

LOG = logging.getLogger(__name__)

tic = time.time()


class SNOfinder:
    """Find Simultaneous Nadir Overpass (SNO) points between two satellites.

    Finding or predicting SNO points between two satellites.
    """

    def __init__(self, platform_id, calipso_id, time_window, sno_min_thr, arc_len_min=2):
        """Initialize the SNO finder class."""
        self.platform_id = platform_id
        self.calipso_id = calipso_id
        self.time_start = time_window[0]
        self.time_end = time_window[1]
        self.arc_len_min = arc_len_min
        self.sno_minute_threshold = sno_min_thr

    def set_configuration(self, configfile):
        """Set the basic configuration from yaml config file."""
        conf = get_config(configfile)
        self._conf = conf

        self.station = {}
        self.station['lon'] = conf['station']['longitude']
        self.station['lat'] = conf['station']['latitude']
        self.station['alt'] = conf['station']['altitude']

    def dataframe2geojson(self, df_):
        """Convert the resulting Pandas dataframe to a Geojson object."""
        gjson = {"type": "FeatureCollection", "features": []}

        # Go through dataframe, append entries to geojson format
        for _, row in df_.iterrows():
            feature = {"type": "Feature", "geometry": {"type": "Point",
                                                       "coordinates": [row['sno_longitude'],
                                                                       row['sno_latitude']]},
                       "properties": {"datetime1": row['satAdatetime'].isoformat(),
                                      "datetime2": row['satBdatetime'].isoformat(),
                                      "tdif_fmin": row['minutes_diff'],
                                      "within_area": row["within_local_reception_area"]}}
            gjson['features'].append(feature)
        self.geojson_results = gjson

    def write_geojson(self, filename):
        """Write the geojson results to file."""
        import json
        with open(filename, 'w') as fp:
            json.dump(self.geojson_results, fp)

    def get_snos_within_time_window(self):
        """Search and retrieve the SNOs inside the time window defined."""
        tle_dirs = self._conf['tle-dirs']
        tle_file_format = self._conf['tle-file-format']

        import sys
        # Check if satellite is supported:
        if OSCAR_NAMES.get(self.platform_id, self.platform_id).upper() not in SATELLITES.keys():
            LOG.error("Platform %s not supported!" % self.platform_id)
            sys.exit()

        pobj = Parser(tle_file_format)
        for tledir in tle_dirs:
            filename_calipso = Path(tledir) / pobj.compose({'platform': self.calipso_id})
            if filename_calipso.exists():
                break

        for tledir in tle_dirs:
            filename_other = Path(tledir) / pobj.compose({'platform': self.platform_id})
            if filename_other.exists():
                break

        TLE_ID_CALIPSO = get_satellite_catalogue_number_from_name(self.calipso_id)
        TLE_ID_OTHER = get_satellite_catalogue_number_from_name(self.platform_id)

        minthr_step = 20  # min less than half an orbit probably
        dtime = timedelta(seconds=60 * minthr_step * 2.0)
        timestep_double = timedelta(seconds=60 * minthr_step * 2.0)
        # make sure the two sat pass the SNO in the same step. We need and overlap of at least half minthr minutes.
        timestep_plus_30s = timedelta(seconds=60 * minthr_step * 1.0 + (self.sno_minute_threshold*0.5)*60 + 30)

        calipso_obj = dt.datetime(1970, 1, 1).replace(tzinfo=timezone.utc)
        other_obj = dt.datetime(1970, 1, 1).replace(tzinfo=timezone.utc)

        tobj = self.time_start
        tle_calipso = None
        tle_the_other_one = None
        tobj_tmp = self.time_start
        i = 0
        t_diff = timedelta(days=1)
        results = []
        while tobj < self.time_end:
            i = i + 1
            if i == 100:
                message = time.time() - tic, "seconds", dtime
                print(message)

            if not tle_calipso or calipso_obj - tobj > t_diff or tobj - calipso_obj > t_diff:
                tle_calipso = get_tle_archive(tobj, str(filename_calipso), TLE_ID_CALIPSO, TLE_BUFFER_CALIPSO)

            if tle_calipso:
                calipso_obj = get_datetime_from_tle(tle_calipso)

            # if tle_calipso:
            #     ts = (tle_calipso.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
            #     calipso_obj = dt.datetime.utcfromtimestamp(ts)
            if not tle_the_other_one or other_obj - tobj > t_diff or tobj - other_obj > t_diff:
                tle_the_other_one = get_tle_archive(tobj, filename_other, TLE_ID_OTHER, TLE_BUFFER_OTHER)
            # if tle_the_other_one:
            #     ts = (tle_the_other_one.epoch - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
            #     other_obj = dt.datetime.utcfromtimestamp(ts)
            if tle_the_other_one:
                other_obj = get_datetime_from_tle(tle_the_other_one)

            calipso = Orbital(self.calipso_id,
                              line1=tle_calipso.line1,
                              line2=tle_calipso.line2)
            the_other_one = Orbital(self.platform_id,
                                    line1=tle_the_other_one.line1,
                                    line2=tle_the_other_one.line2)

            got_intersection_acurate = False
            arc_calipso_vector = get_arc_vector(tobj, timestep_plus_30s, calipso, self.arc_len_min)
            arc_the_other_one_vector = get_arc_vector(tobj, timestep_plus_30s, the_other_one,  self.arc_len_min)
            # Approximate tracks with one arc each self.arc_len_min minutes.
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
                                    tobj, self.sno_minute_threshold, self.station)

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
                LOG.debug(tobj_tmp.strftime("%Y-%m-%d"))

        print(str(results[0]))
        return pd.DataFrame(results)


def get_satellite_catalogue_number_from_name(satname):
    """Return the satellite catalogue number from the platform name."""
    return SATELLITES[OSCAR_NAMES.get(satname, satname).upper()]


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


def get_sno_point(calipso, the_other_one, arc_calipso, arc_the_other_one, tobj, minthr, station):
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
                                                station['lon'], station['lat'], station['alt'])
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


def get_closest_sno_to_reference(all_features, rfeature, tol_seconds=ZERO_SECONDS):
    """Check all features found and find the one matching the reference.

    Returns the distance between the two SNOs in km.
    """
    rgeom = rfeature['geometry']
    rlon, rlat = rgeom['coordinates']
    rprop = rfeature['properties']

    rtime1 = dt.datetime.fromisoformat(rprop['datetime1'])
    rtime2 = dt.datetime.fromisoformat(rprop['datetime2'])

    overlap_found = False
    for tfeat in all_features['features']:
        tgeom = tfeat['geometry']
        tlon, tlat = tgeom['coordinates']
        tprop = tfeat['properties']

        ttime1 = dt.datetime.fromisoformat(tprop['datetime1'])
        ttime2 = dt.datetime.fromisoformat(tprop['datetime2'])

        # Check for time overlap:
        if check_overlapping_times((rtime1, rtime2), (ttime1, ttime2), tol_seconds):
            # print("Overlap in times")
            if tol_seconds > ZERO_SECONDS:
                print("Approximate overlap in times")
            km_dist = distance.distance((tlat, tlon), (rlat, rlon)).kilometers
            overlap_found = True
            print(f"Distance between SNOs: {km_dist:5.1f} km")
            break

    if overlap_found:
        return km_dist, (rtime1, rtime2), (ttime1, ttime2)
    else:
        return None, (rtime1, rtime2), None


def geojson_compare_derived_snos_against_reference(this_features, ref_features):
    """Compare the SNOs derived against a reference set.

    We loop over the Geojson features of the reference dataset and see if a
    corresponding SNO can be found among the Geojson features output of the SNO
    finder.
    """
    for rfeat in ref_features['features']:
        km_dist, twindow_ref, twindow_this = get_closest_sno_to_reference(this_features, rfeat)
        if not km_dist:
            rgeom = rfeat['geometry']
            rlon, rlat = rgeom['coordinates']

            print("Failed to find strict overlap in times")
            km_dist, twindow_ref, twindow_this = get_closest_sno_to_reference(this_features, rfeat,
                                                                              timedelta(seconds=10))
            if not km_dist:
                print("Failed to find overlap in times")


def check_overlapping_times(twindow1, twindow2, tol_seconds=ZERO_SECONDS):
    """Check if two time windows overlap.

    A tolerance can be given to allow accepting time windows that are very
    close but not actuall overlapping.
    """
    twindow1 = check_time_window(twindow1)
    twindow2 = check_time_window(twindow2)

    rtime1 = twindow1[0] - tol_seconds
    rtime2 = twindow1[1] + tol_seconds
    ttime1 = twindow2[0] - tol_seconds
    ttime2 = twindow2[1] + tol_seconds

    if (rtime1 >= ttime1 and rtime1 <= ttime2) or \
       (rtime2 >= ttime1 and rtime2 <= ttime2):
        return True
    return False


def check_time_window(twindow):
    """Check the consistency of a time window, so that the first item is always the older."""
    time1 = twindow[0]
    time2 = twindow[1]
    if time1 > time2:
        tmp_time = time1
        time1 = time2
        time2 = tmp_time
        twindow = (time1, time2)

    return twindow
