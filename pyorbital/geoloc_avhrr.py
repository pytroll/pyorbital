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

"""This module provides geoloc operations specific to the avhrr instrument."""

import logging

import numpy as np
from pyproj import Geod

from pyorbital.geoloc import ScanGeometry, compute_pixels, get_lonlatalt

logger = logging.getLogger(__name__)
geod = Geod(ellps="WGS84")

def compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, start_time, tle) -> None:
    """Compute the longitute, latitude and altitude of given gcps (lines, columns of the swath)."""
    # TODO: account for missing scanlines
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
    """Estimate time offset and attitude deviations from gcps."""
    from scipy.optimize import minimize

    # we need to work in seconds*1e3 to avoid the nanosecond precision issue
    res = minimize(compute_gcp_distances_to_reference_lonlats,
                   x0=(0, 0, 0, 0),
                   args=(gcps, start_time, tle, max_scan_angle, (ref_lons, ref_lats)),
                   bounds=((-0.02, 0.02) , (-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)))
    if not res.success:
        raise RuntimeError("Time and attitude estimation did not converge")
    time_diff, roll, pitch, yaw = res.x * [1e3, 1, 1, 1]
    logger.debug(f"Estimated time difference to {time_diff} seconds, attitude to {roll}, {pitch}, {yaw} degrees")

    return time_diff, roll, pitch, yaw



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
    return np.sum(distances**2)
