#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2025 PyTroll Community

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

"""This module provides geoloc operations specific to the avhrr instrument.

In particular, it provides functions that allow matching gcp location in swath coordinates to reference positions, and
then minimise the distance to these positions by adjusting the time offset and attitude.
"""

import logging

import numpy as np
from pyproj import Geod

from pyorbital.geoloc import ScanGeometry, compute_pixels, get_lonlatalt

logger = logging.getLogger(__name__)
geod = Geod(ellps="WGS84")

def compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, start_time, tle) -> None:
    """Compute the longitute, latitude and altitude of given gcps (scanlines, columns of the swath).

    The gcps are arbitrary location in swath coordinates, for example (10.3, 7.7) for a gcp at line 10.3 in the swath,
    and column 7.7. This function returns the geographical coordinates of the gcps.

    The scanlines are relative to the pass scanline numbers, zero-based.
    """
    time_line_interval = 1/6
    time_row_interval = 25e-6

    fov_x = gcps[:, 1]
    fov_y = gcps[:, 0]

    scan_angles_across = (fov_x / 1023.5 - 1) * np.deg2rad(-max_scan_angle)
    scan_angles_along = np.zeros_like(scan_angles_across)
    scan_angles = np.vstack((scan_angles_across, scan_angles_along))
    time_offsets = np.array(fov_x * time_row_interval + fov_y * time_line_interval)
    geom = ScanGeometry(scan_angles, time_offsets)
    start_time = np.datetime64(start_time)
    s_times = geom.times(start_time)

    pixels_pos = compute_pixels(tle, geom, s_times, rpy)
    return get_lonlatalt(pixels_pos, s_times)


def estimate_time_and_attitude_deviations(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle):
    """Estimate time offset and attitude deviations from gcps.

    Provided reference longitudes and latitudes for the gcps, this function minimises the attitude and time offset
    needed to match the gcp coordinates to the reference coordinates.
    """
    from scipy.optimize import minimize

    distances = compute_gcp_distances_to_reference_lonlats((0, 0, 0, 0), gcps, start_time, tle, max_scan_angle,
                                                           (ref_lons, ref_lats))
    logger.debug(f"GCP distances: median {np.median(distances)}, std {np.std(distances)}")
    # we need to work in seconds*1e3 to avoid the nanosecond precision issue
    res = minimize(compute_gcp_accumulated_squared_distances_to_reference_lonlats,
                   x0=(0, 0, 0, 0),
                   args=(gcps, start_time, tle, max_scan_angle, (ref_lons, ref_lats)),
                   bounds=((-0.007, 0.007) , (-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)))
    if not res.success:
        raise RuntimeError("Time and attitude estimation did not converge")
    time_diff, roll, pitch, yaw = res.x * [1e3, 1, 1, 1]
    logger.debug(f"Estimated time difference to {time_diff} seconds, attitude to {roll}, {pitch}, {yaw} degrees")
    distances = compute_gcp_distances_to_reference_lonlats(res.x, gcps, start_time, tle, max_scan_angle,
                                                           (ref_lons, ref_lats))
    logger.debug(f"Remaining GCP distances: median {np.median(distances)}, std {np.std(distances)}")

    return time_diff, roll, pitch, yaw


def compute_gcp_accumulated_squared_distances_to_reference_lonlats(
        variables, gcps, start_time, tle, max_scan_angle, refs):
    """Compute the summed squared distance fot gcps to reference lonlats.

    Given the gcps (in swath coordinates) along with attitude and time offset, compute the sum of squared distances to
    the reference lons and lats of the gcps.
    """
    distances = compute_gcp_distances_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs)
    return np.sum(distances**2)


def compute_gcp_distances_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs):
    """Compute the gcp distances to references lonlats."""
    time_diff, roll, pitch, yaw = variables
    # we need to work in seconds*1e3 to avoid the nanosecond precision issue
    time = np.datetime64(start_time) + np.timedelta64(int(time_diff * 1e12), "ns")
    lons, lats, _ = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, (roll, pitch, yaw), time, tle)
    valid = np.isfinite(lons)
    lons = lons[valid]
    lats = lats[valid]
    ref_lons, ref_lats = refs
    ref_lons = np.array(ref_lons)[valid]
    ref_lats = np.array(ref_lats)[valid]
    _, _, distances = geod.inv(ref_lons, ref_lats, lons, lats)
    return distances
