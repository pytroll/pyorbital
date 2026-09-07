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

def compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, start_time, tle, yaw_steering=False,
                                 nadir_convention=None) -> None:
    """Compute the longitute, latitude and altitude of given gcps (scanlines, columns of the swath).

    The gcps are arbitrary location in swath coordinates, for example (10.3, 7.7) for a gcp at line 10.3 in the swath,
    and column 7.7. This function returns the geographical coordinates of the gcps.

    The scanlines are relative to the pass scanline numbers, zero-based.

    Pass *yaw_steering* for a platform that turns as it flies to hold its swath
    square to the ground track, as Metop does and the POES platforms do not. It
    must match the convention the geolocation under study was computed with,
    since a mismatch shows up as a whole-swath yaw of a few degrees.

    *nadir_convention* must likewise match, for the same reason: a model standing
    on a different nadir than the navigation it is fitted to absorbs the
    difference, which reaches some hundreds of metres at the swath edge.
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

    pixels_pos = compute_pixels(tle, geom, s_times, rpy, yaw_steering=yaw_steering,
                                nadir_convention=nadir_convention)
    return get_lonlatalt(pixels_pos, s_times)


# The minimiser carries the time offset in kiloseconds so that a step it considers
# small is still far larger than the nanosecond the timestamps are stored in.
TIME_SEARCH_REACH = 0.007   # kiloseconds, so seven seconds either side
ATTITUDE_SEARCH_REACH = 0.5  # radians on each angle, about 28 degrees


def _with_time(searched, solve_for_time):
    """Return the full four variables, putting back the time when it was not searched for."""
    return np.asarray(searched) if solve_for_time else np.concatenate([[0.0], searched])


def estimate_time_and_attitude_deviations(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle,
                                          yaw_steering=False, nadir_convention=None,
                                          time_offset_guess=0.0, solve_for_time=True):
    """Estimate time offset and attitude deviations from gcps.

    Provided reference longitudes and latitudes for the gcps, this function minimises the attitude and time offset
    needed to match the gcp coordinates to the reference coordinates.

    The search reaches only seven seconds, which is a deliberate guard: a time offset the
    data cannot pin down would otherwise wander off and drag the attitude with it, since
    a shift along the track can be written either as time or as pitch. When something
    upstream already knows roughly how far the swath has moved -- a coarse image match,
    say -- pass that as *time_offset_guess* in seconds, and the search reaches seven
    seconds either side of it rather than either side of zero.

    Platforms whose clock is disciplined -- the KLM series and Metop, whose along-track
    displacement never leaves the coarse matcher's noise floor -- should pass
    *solve_for_time* as False. Solving for an offset already known to be zero only
    lets the pitch absorb its noise, since the two are barely distinguishable.
    """
    from scipy.optimize import least_squares, minimize

    original_distances = compute_gcp_distances_to_reference_lonlats(
        (0, 0, 0, 0), gcps, start_time, tle, max_scan_angle, (ref_lons, ref_lats),
        yaw_steering, nadir_convention)
    original_median_distance = np.median(original_distances)
    logger.debug(f"GCP distances: median {original_median_distance}, std {np.std(original_distances)}")
    guessed = time_offset_guess / 1e3
    reach = np.array((TIME_SEARCH_REACH, ATTITUDE_SEARCH_REACH, ATTITUDE_SEARCH_REACH,
                      ATTITUDE_SEARCH_REACH))
    middle = np.array((guessed, 0.0, 0.0, 0.0))
    held = slice(None) if solve_for_time else slice(1, None)

    def offsets(searched, *args):
        return compute_gcp_offsets_to_reference_lonlats(_with_time(searched, solve_for_time), *args)

    res = least_squares(offsets,
                        x0=middle[held],
                        args=(gcps, start_time, tle, max_scan_angle, (ref_lons, ref_lats), yaw_steering,
                              nadir_convention),
                        bounds=(middle[held] - reach[held], middle[held] + reach[held]), x_scale="jac")
    if not res.success:
        raise RuntimeError("Time and attitude estimation did not converge")
    settled = _with_time(res.x, solve_for_time)
    if solve_for_time and res.active_mask[0] != 0:
        raise RuntimeError("The time offset did not settle inside its search; "
                           "nothing in the data holds it, and the attitude pays for it")
    time_diff, roll, pitch, yaw = settled * [1e3, 1, 1, 1]
    logger.debug(f"Estimated time difference to {time_diff} seconds, "
                 f"attitude to {np.rad2deg(roll)}, {np.rad2deg(pitch)}, {np.rad2deg(yaw)} degrees")
    distances = compute_gcp_distances_to_reference_lonlats(settled, gcps, start_time, tle, max_scan_angle,
                                                           (ref_lons, ref_lats), yaw_steering, nadir_convention)

    minimized_median_distance = np.median(distances)
    logger.debug(f"Remaining GCP distances: median {minimized_median_distance}, std {np.std(distances)}")

    return time_diff, (roll, pitch, yaw), (original_distances, distances)


def estimate_time_offset(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle):
    """Estimate time offset from gcps.

    Provided reference longitudes and latitudes for the gcps, this function minimises the time offset
    needed to match the gcp coordinates to the reference coordinates.
    """
    from scipy.optimize import minimize

    original_distances = compute_gcp_distances_to_reference_lonlats((0, 0, 0, 0), gcps, start_time, tle, max_scan_angle,
                                                                    (ref_lons, ref_lats))
    original_median_distance = np.median(original_distances)
    logger.debug(f"GCP distances: median {original_median_distance}, std {np.std(original_distances)}")

    def gcp_distance_for_time(time):
        dist = compute_gcp_accumulated_squared_distances_to_reference_lonlats((time[0], 0, 0, 0), gcps, start_time, tle,
                                                                              max_scan_angle, (ref_lons, ref_lats))
        return dist

    # we need to work in seconds*1e3 to avoid the nanosecond precision issue
    res = minimize(gcp_distance_for_time,
                   x0=(0,),
                   bounds=((-0.03, 0.03),),
                   options=dict(ftol=1e-1),
                   )
    if not res.success:
        raise RuntimeError("Time offset estimation did not converge")
    time_diff, = res.x * [1e3,]
    logger.debug(f"Estimated time difference to {time_diff} seconds")
    distances = compute_gcp_distances_to_reference_lonlats((res.x[0], 0, 0, 0), gcps, start_time, tle, max_scan_angle,
                                                           (ref_lons, ref_lats))

    minimized_median_distance = np.median(distances)
    logger.debug(f"Remaining GCP distances: median {minimized_median_distance}, std {np.std(distances)}")

    return time_diff, (original_distances, distances)


def compute_gcp_accumulated_squared_distances_to_reference_lonlats(
        variables, gcps, start_time, tle, max_scan_angle, refs, yaw_steering=False,
        nadir_convention=None):
    """Compute the summed squared distance fot gcps to reference lonlats.

    Given the gcps (in swath coordinates) along with attitude and time offset, compute the sum of squared distances to
    the reference lons and lats of the gcps.
    """
    distances = compute_gcp_distances_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs,
                                                           yaw_steering, nadir_convention)
    return np.sum(distances**2)


def _misses_from_reference(variables, gcps, start_time, tle, max_scan_angle, refs,
                           yaw_steering=False, nadir_convention=None):
    """Return which way and how far each gcp landed from its reference point."""
    time_diff, roll, pitch, yaw = variables
    time = np.datetime64(start_time) + np.timedelta64(int(time_diff * 1e12), "ns")
    lons, lats, _ = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, (roll, pitch, yaw), time, tle,
                                                yaw_steering, nadir_convention)
    valid = np.isfinite(lons)
    lons = lons[valid]
    lats = lats[valid]
    ref_lons, ref_lats = refs
    ref_lons = np.array(ref_lons)[valid]
    ref_lats = np.array(ref_lats)[valid]
    bearings, _, distances = geod.inv(ref_lons, ref_lats, lons, lats)
    return np.radians(bearings), distances


def compute_gcp_distances_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs,
                                               yaw_steering=False, nadir_convention=None):
    """Compute the gcp distances to references lonlats."""
    _, distances = _misses_from_reference(variables, gcps, start_time, tle, max_scan_angle, refs,
                                          yaw_steering, nadir_convention)
    return distances


def compute_gcp_offsets_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs,
                                             yaw_steering=False, nadir_convention=None):
    """Return each gcp's miss as a northward and an eastward component, in metres.

    The same misses as the distances, kept signed instead of collapsed into a magnitude.
    A least-squares solver wants the residuals themselves; a magnitude has already thrown
    away which way each point missed, which is most of what tells the parameters apart.
    """
    bearings, distances = _misses_from_reference(variables, gcps, start_time, tle, max_scan_angle, refs,
                                                 yaw_steering, nadir_convention)
    return np.concatenate([distances * np.cos(bearings), distances * np.sin(bearings)])
