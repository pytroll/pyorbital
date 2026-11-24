"""This module provides geoloc operations specific to the avhrr instrument.

In particular, it provides functions that allow matching gcp location in swath coordinates to reference positions, and
then minimise the distance to these positions by adjusting the time offset and attitude.
"""

import logging
import warnings
from typing import NamedTuple

import numpy as np
from pyproj import Geod

from pyorbital.geoloc import ScanGeometry, _resolve_convention, compute_pixels, get_lonlatalt

logger = logging.getLogger(__name__)
geod = Geod(ellps="WGS84")

def compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, start_time, tle, yaw_steering=False,
                                 nadir_convention=None, rotation_order=None) -> None:
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
                                nadir_convention=nadir_convention,
                                rotation_order=rotation_order)
    return get_lonlatalt(pixels_pos, s_times)


# The minimiser carries the time offset in kiloseconds so that a step it considers
# small is still far larger than the nanosecond the timestamps are stored in.
TIME_SEARCH_REACH = 0.007   # kiloseconds, so seven seconds either side
ATTITUDE_SEARCH_REACH = 0.5  # radians on each angle, about 28 degrees


def _all_four(searched, free):
    """Return the full four variables, holding at zero the ones not searched for.

    A parameter that was not searched for is not unknown: it is being asserted to be
    zero, which is the whole point of holding it.
    """
    full = np.zeros(4)
    full[free] = searched
    return full


def _search_box(time_offset_guess, solve_for_time, solve_for_pitch):
    """Return the centre and half-width of the search, and which of the four parameters it frees."""
    reach = np.array((TIME_SEARCH_REACH, ATTITUDE_SEARCH_REACH, ATTITUDE_SEARCH_REACH,
                      ATTITUDE_SEARCH_REACH))
    middle = np.array((time_offset_guess / 1e3, 0.0, 0.0, 0.0))
    free = np.array([solve_for_time, True, solve_for_pitch, True])
    return middle, reach, free


def _refuse_an_unsettled_time(res, solve_for_time):
    """Refuse a time offset that came to rest on the wall of its search."""
    if solve_for_time and res.active_mask[0] != 0:
        raise RuntimeError("The time offset did not settle inside its search; "
                           "nothing in the data holds it, and the attitude pays for it")


def estimate_time_and_attitude_deviations(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle,
                                          yaw_steering=False, nadir_convention=None,
                                          rotation_order=None, solve_for_pitch=True,
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

    A shift along the track can be written as a time offset or as a pitch, and the two
    exchange at about 2.25 seconds per degree. Only the curvature the shift leaves
    across the swath tells them apart, and a single pass rarely constrains that, so
    solving for both together lets each absorb the other's noise. Pass
    *solve_for_time* or *solve_for_pitch* as False to hold one of them at zero and
    fit the other. Which to hold is a property of the platform: one whose clock is
    steered in flight can be believed and its pitch fitted, while one whose clock
    drifts must have its time fitted instead.

    A parameter held this way is not being ignored. It is being asserted to be zero,
    and the answer returned for it is that assertion rather than a measurement.

    .. deprecated::
        Use :func:`fit_navigation`, which judges each miss in swath pixels rather than
        metres, keeps a few bad matches from dragging the rest, and returns each
        point's leftover miss and the uncertainty of what it found.
    """
    warnings.warn(
        "estimate_time_and_attitude_deviations is deprecated and will be removed; use "
        "fit_navigation, which judges each miss in swath pixels, robustly, and says how "
        "far its answer can be trusted",
        DeprecationWarning, stacklevel=2)
    from scipy.optimize import least_squares

    original_distances = compute_gcp_distances_to_reference_lonlats(
        (0, 0, 0, 0), gcps, start_time, tle, max_scan_angle, (ref_lons, ref_lats),
        yaw_steering, nadir_convention, rotation_order)
    original_median_distance = np.median(original_distances)
    logger.debug(f"GCP distances: median {original_median_distance}, std {np.std(original_distances)}")
    middle, reach, free = _search_box(time_offset_guess, solve_for_time, solve_for_pitch)

    def offsets(searched, *args):
        return compute_gcp_offsets_to_reference_lonlats(_all_four(searched, free), *args)

    res = least_squares(offsets,
                        x0=middle[free],
                        args=(gcps, start_time, tle, max_scan_angle, (ref_lons, ref_lats), yaw_steering,
                              nadir_convention, rotation_order),
                        bounds=(middle[free] - reach[free], middle[free] + reach[free]), x_scale="jac")
    if not res.success:
        raise RuntimeError("Time and attitude estimation did not converge")
    settled = _all_four(res.x, free)
    _refuse_an_unsettled_time(res, solve_for_time)
    time_diff, roll, pitch, yaw = settled * [1e3, 1, 1, 1]
    logger.debug(f"Estimated time difference to {time_diff} seconds, "
                 f"attitude to {np.rad2deg(roll)}, {np.rad2deg(pitch)}, {np.rad2deg(yaw)} degrees")
    distances = compute_gcp_distances_to_reference_lonlats(settled, gcps, start_time, tle, max_scan_angle,
                                                           (ref_lons, ref_lats), yaw_steering, nadir_convention,
                                                           rotation_order)

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
        nadir_convention=None, rotation_order=None):
    """Compute the summed squared distance fot gcps to reference lonlats.

    Given the gcps (in swath coordinates) along with attitude and time offset, compute the sum of squared distances to
    the reference lons and lats of the gcps.
    """
    distances = compute_gcp_distances_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs,
                                                           yaw_steering, nadir_convention, rotation_order)
    return np.sum(distances**2)


def _misses_from_reference(variables, gcps, start_time, tle, max_scan_angle, refs,
                           yaw_steering=False, nadir_convention=None, rotation_order=None):
    """Return which way and how far each gcp landed from its reference point."""
    time_diff, roll, pitch, yaw = variables
    time = np.datetime64(start_time) + np.timedelta64(int(time_diff * 1e12), "ns")
    lons, lats, _ = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, (roll, pitch, yaw), time, tle,
                                                yaw_steering, nadir_convention, rotation_order)
    valid = np.isfinite(lons)
    lons = lons[valid]
    lats = lats[valid]
    ref_lons, ref_lats = refs
    ref_lons = np.array(ref_lons)[valid]
    ref_lats = np.array(ref_lats)[valid]
    bearings, _, distances = geod.inv(ref_lons, ref_lats, lons, lats)
    return np.radians(bearings), distances


def compute_gcp_distances_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs,
                                               yaw_steering=False, nadir_convention=None,
                                               rotation_order=None):
    """Compute the gcp distances to references lonlats."""
    _, distances = _misses_from_reference(variables, gcps, start_time, tle, max_scan_angle, refs,
                                          yaw_steering, nadir_convention, rotation_order)
    return distances


def compute_gcp_offsets_to_reference_lonlats(variables, gcps, start_time, tle, max_scan_angle, refs,
                                             yaw_steering=False, nadir_convention=None,
                                             rotation_order=None):
    """Return each gcp's miss as a northward and an eastward component, in metres.

    The same misses as the distances, kept signed instead of collapsed into a magnitude.
    A least-squares solver wants the residuals themselves; a magnitude has already thrown
    away which way each point missed, which is most of what tells the parameters apart.
    """
    bearings, distances = _misses_from_reference(variables, gcps, start_time, tle, max_scan_angle, refs,
                                                 yaw_steering, nadir_convention, rotation_order)
    return np.concatenate([distances * np.cos(bearings), distances * np.sin(bearings)])


def _north_east(from_lons, from_lats, to_lons, to_lats):
    """Return the step from one set of points to another, as northward and eastward metres."""
    bearings, _, distances = geod.inv(from_lons, from_lats, to_lons, to_lats)
    bearings = np.radians(bearings)
    return np.stack([distances * np.cos(bearings), distances * np.sin(bearings)], axis=-1)


def pixel_metric(variables, gcps, start_time, tle, max_scan_angle, yaw_steering=False,
                 nadir_convention=None, rotation_order=None):
    """Return, for each gcp, the matrix turning a ground step in metres into swath pixels.

    A step of one line and one sample is taken from every gcp and its length and direction
    on the ground measured; inverting those two steps gives the conversion. The result is
    a stack of 2x2 matrices taking (north, east) metres to (lines, samples).

    It holds only near the geometry it was computed for, but the corrections a fit makes
    are small fractions of a degree, so the footprint barely changes across them and the
    conversion can be computed once and kept.
    """
    time_diff, roll, pitch, yaw = variables
    time = np.datetime64(start_time) + np.timedelta64(int(time_diff * 1e12), "ns")
    gcps = np.asarray(gcps, dtype=float)
    stepped = np.concatenate([gcps, gcps + [1.0, 0.0], gcps + [0.0, 1.0]])
    lons, lats, _ = compute_avhrr_gcps_lonlatalt(stepped, max_scan_angle, (roll, pitch, yaw), time, tle,
                                                yaw_steering, nadir_convention, rotation_order)
    lons, lats = np.split(lons, 3), np.split(lats, 3)
    per_line = _north_east(lons[0], lats[0], lons[1], lats[1])
    per_sample = _north_east(lons[0], lats[0], lons[2], lats[2])
    metres_per_pixel = np.stack([per_line, per_sample], axis=-1)
    usable = np.all(np.isfinite(metres_per_pixel), axis=(1, 2))
    metric = np.full_like(metres_per_pixel, np.nan)
    metric[usable] = np.linalg.inv(metres_per_pixel[usable])
    return metric


def compute_gcp_misses_in_pixels(variables, gcps, start_time, tle, max_scan_angle, refs,
                                 yaw_steering=False, nadir_convention=None, rotation_order=None,
                                 metric=None):
    """Return how far each gcp lands from its reference, in swath lines and samples.

    The miss is where the navigation puts the point less where it truly is, so a
    reference one sample further along the scan than the gcp gives a miss of minus one
    sample. A gcp that cannot be placed keeps its row, as NaN, so the misses stay
    aligned with the gcps they belong to.

    Pass a *metric* from :func:`pixel_metric` to reuse one conversion across many calls;
    without one it is computed at *variables*.
    """
    if metric is None:
        metric = pixel_metric(variables, gcps, start_time, tle, max_scan_angle, yaw_steering,
                              nadir_convention, rotation_order)
    time_diff, roll, pitch, yaw = variables
    time = np.datetime64(start_time) + np.timedelta64(int(time_diff * 1e12), "ns")
    lons, lats, _ = compute_avhrr_gcps_lonlatalt(np.asarray(gcps, dtype=float), max_scan_angle,
                                                (roll, pitch, yaw), time, tle,
                                                yaw_steering, nadir_convention, rotation_order)
    ref_lons, ref_lats = (np.asarray(ref, dtype=float) for ref in refs)
    misses = _north_east(ref_lons, ref_lats, lons, lats)
    return np.einsum("nij,nj->ni", metric, misses)


class NavigationFit(NamedTuple):
    """What a navigation fit found, and how well it accounts for each control point.

    *time_offset* is in seconds and *attitude* is roll, pitch and yaw in radians.
    *misses* are what is left of each gcp's miss after the fit, in swath lines and
    samples, one row per gcp in the order given. *distances* are the same leftovers
    measured on the ground, in metres, for products that already quote them so.

    *covariance* is that of time offset, roll, pitch and yaw, in seconds and radians.
    A parameter that was held rather than fitted has no row or column in it beyond
    zeros: it was asserted, not measured.

    *settings* records every choice the fit was made with, conventions resolved to
    the values actually used, so that a later fit can reproduce or extend it.
    """

    time_offset: float
    attitude: tuple
    misses: np.ndarray
    distances: np.ndarray
    covariance: np.ndarray
    settings: dict

    @property
    def uncertainty(self):
        """Return the standard deviation of time offset, roll, pitch and yaw, in seconds and radians."""
        return np.sqrt(np.diag(self.covariance))

    @property
    def correlation(self):
        """Return how the four parameters' errors move together, NaN where one was held.

        Time and pitch are nearly the same observable along the track, so their
        correlation says whether the split between them can be believed.
        """
        spread = self.uncertainty
        with np.errstate(divide="ignore", invalid="ignore"):
            correlation = self.covariance / np.outer(spread, spread)
        correlation[:, spread == 0] = np.nan
        correlation[spread == 0, :] = np.nan
        return correlation


def fit_navigation(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle,
                   yaw_steering=False, nadir_convention=None, rotation_order=None,
                   solve_for_pitch=True, time_offset_guess=0.0, solve_for_time=True,
                   loss="soft_l1", loss_scale=1.0, heights=None):
    """Fit the time offset and attitude that carry the gcps onto their references.

    Each miss is judged in swath pixels rather than metres. The matching that found the
    gcps is about equally precise, in pixels, anywhere in the swath, whereas a pixel at
    the swath edge covers several times the ground it does at nadir: fitted in metres,
    the edge would outweigh the centre for no better reason than its footprint.

    Squared misses would let each gcp pull in proportion to how far out it is, so a
    few confident mismatches would drag every good point towards them. The *loss* is
    therefore robust by default: misses within *loss_scale* pixels count as squares,
    those beyond it roughly in proportion, which caps how hard any one point can pull.
    Any loss :func:`scipy.optimize.least_squares` takes may be given, ``"linear"``
    for plain least squares.

    A reference stands at its true place, but the swath shows it where the line of sight
    to it meets the ellipsoid, which for anything above sea level is further from nadir.
    Given the *heights* of the references in metres, each is moved to where it is seen
    before the fit, so the terrain is not taken for navigation error and the fitted
    geometry stays that of the ellipsoid. Without them, the references are taken to lie
    on the ellipsoid; a missing height counts as sea level.

    The search reach and the parameters that may be held are as for
    :func:`estimate_time_and_attitude_deviations`.
    """
    from scipy.optimize import least_squares

    nadir_convention = _resolve_convention(nadir_convention, "nadir_convention", "nadir convention",
                                           "nadir_convention='geocentric'")
    rotation_order = _resolve_convention(rotation_order, "rotation_order", "rotation order",
                                         "rotation_order='pitch_first'")
    settings = {"solve_for_time": bool(solve_for_time), "solve_for_pitch": bool(solve_for_pitch),
                "time_offset_guess": float(time_offset_guess), "yaw_steering": bool(yaw_steering),
                "nadir_convention": nadir_convention, "rotation_order": rotation_order,
                "max_scan_angle": float(max_scan_angle), "loss": loss, "loss_scale": float(loss_scale),
                "residual_units": "swath pixels (lines, samples)"}
    if heights is not None:
        ref_lons, ref_lats = where_the_references_are_seen(gcps, ref_lons, ref_lats, heights, start_time, tle,
                                                           time_offset_guess)
    settings["terrain_heights"] = heights is not None
    refs = (ref_lons, ref_lats)
    geometry = (yaw_steering, nadir_convention, rotation_order)
    middle, reach, free = _search_box(time_offset_guess, solve_for_time, solve_for_pitch)
    metric = pixel_metric(middle, gcps, start_time, tle, max_scan_angle, *geometry)

    def misses(variables):
        return compute_gcp_misses_in_pixels(variables, gcps, start_time, tle, max_scan_angle, refs,
                                            *geometry, metric=metric)

    def residuals(searched):
        flat = misses(_all_four(searched, free)).ravel()
        return flat[np.isfinite(flat)]

    res = least_squares(residuals, x0=middle[free],
                        bounds=(middle[free] - reach[free], middle[free] + reach[free]), x_scale="jac",
                        loss=loss, f_scale=loss_scale)
    if not res.success:
        raise RuntimeError("Time and attitude estimation did not converge")
    _refuse_an_unsettled_time(res, solve_for_time)
    settled = _all_four(res.x, free)
    time_diff, roll, pitch, yaw = settled * [1e3, 1, 1, 1]
    left = misses(settled)
    placed = left[np.all(np.isfinite(left), axis=1)]
    lines_rms, samples_rms = np.sqrt(np.mean(placed ** 2, axis=0)) if len(placed) else (np.nan, np.nan)
    judged = judge_what_a_fit_left(gcps, left)
    logger.debug(f"Fitted {len(placed)} control points: time offset {time_diff:.2f} s, "
                 f"attitude {np.rad2deg(roll):.4f}, {np.rad2deg(pitch):.4f}, {np.rad2deg(yaw):.4f} degrees "
                 f"(roll, pitch, yaw); left over, rms {lines_rms:.2f} lines, {samples_rms:.2f} samples; "
                 f"{judged.mismatches} mismatches; "
                 f"trusted rms {judged.trusted_rms[0]:.2f} lines, {judged.trusted_rms[1]:.2f} samples; "
                 f"smooth part {judged.systematic_rms[0]:.2f} lines, {judged.systematic_rms[1]:.2f} samples")
    return NavigationFit(time_diff, (roll, pitch, yaw), left,
                         _ground_distances(settled, gcps, start_time, tle, max_scan_angle, refs, *geometry),
                         _covariance(res, free), settings)


def _ground_distances(variables, gcps, start_time, tle, max_scan_angle, refs, *geometry):
    """Return each gcp's distance from its reference in metres, NaN in place of any not placed."""
    time_diff, roll, pitch, yaw = variables
    time = np.datetime64(start_time) + np.timedelta64(int(time_diff * 1e12), "ns")
    lons, lats, _ = compute_avhrr_gcps_lonlatalt(np.asarray(gcps, dtype=float), max_scan_angle,
                                                (roll, pitch, yaw), time, tle, *geometry)
    return np.hypot(*_north_east(*(np.asarray(ref, dtype=float) for ref in refs), lons, lats).T)


#: Scales a median absolute deviation to the standard deviation it stands for under normal noise.
_MAD_TO_SIGMA = 1.4826


def _covariance(res, free):
    """Return the covariance of the four parameters from a finished least-squares fit.

    The spread of the misses is read from their median absolute deviation rather than
    their sum of squares, so a few bad matches the loss already discounted do not
    inflate it again. Time is converted from the kiloseconds the solver searched in.
    """
    residuals = res.fun
    dof = max(residuals.size - res.x.size, 1)
    scatter = _MAD_TO_SIGMA * np.median(np.abs(residuals - np.median(residuals)))
    variance = scatter ** 2 * residuals.size / dof
    searched = variance * np.linalg.pinv(res.jac.T @ res.jac)
    covariance = np.zeros((4, 4))
    covariance[np.ix_(free, free)] = searched
    units = np.array([1e3, 1.0, 1.0, 1.0])
    return covariance * np.outer(units, units)


#: How far, in robust standard deviations of its two axes taken together, a miss may stray
#: from the rest before it is taken for a mismatch rather than a measurement. Under normal
#: matching noise one good point in a thousand strays further (the 99.9 % radius of two
#: independent axes); judging each axis alone at three would set aside five in a thousand.
MISMATCH_THRESHOLD = 3.72


class WhatAFitLeft(NamedTuple):
    """How far the trustworthy control points ended up, and how much of it is navigation.

    *trusted_rms* is the root mean square of the misses left on the points that are not
    mismatches, in swath lines and samples. *mismatches* counts the points set aside.
    *systematic_rms* is the part of those misses that varies smoothly across the swath,
    in lines and samples: navigation the fit did not remove, as opposed to the matcher's
    own noise, which differs from point to point.
    """

    trusted_rms: tuple
    mismatches: int
    systematic_rms: tuple


def judge_what_a_fit_left(gcps, misses, threshold=MISMATCH_THRESHOLD):
    """Judge the misses a navigation fit left on its control points.

    A point is a mismatch when its miss strays more than *threshold* from the median of
    all the misses, each axis measured in its own robust standard deviation, read from
    the median absolute deviation. The rest are trusted.

    The smooth part is fitted to the trusted misses with the terms a smooth geometry
    leaves -- a constant, a slope along and across the swath, and their product -- and
    quoted as its root mean square over those points. It cannot fully tell noise from
    geometry: pure noise of spread s leaves about s times the square root of four over
    the number of points in it.
    """
    gcps = np.asarray(gcps, dtype=float)
    misses = np.asarray(misses, dtype=float)
    placed = np.all(np.isfinite(misses), axis=1) & np.all(np.isfinite(gcps), axis=1)
    gcps, misses = gcps[placed], misses[placed]
    if len(misses) == 0:
        return WhatAFitLeft((np.nan, np.nan), 0, (np.nan, np.nan))

    centre = np.median(misses, axis=0)
    spread = _MAD_TO_SIGMA * np.median(np.abs(misses - centre), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        strays = np.hypot(*np.where(spread > 0, (misses - centre) / spread, 0.0).T)
    trusted = strays <= threshold
    kept, where = misses[trusted], gcps[trusted]
    trusted_rms = tuple(float(value) for value in np.sqrt(np.mean(kept ** 2, axis=0)))

    line, sample = ((where - where.mean(axis=0)) / np.maximum(where.std(axis=0), 1.0)).T
    smooth = np.column_stack([np.ones_like(line), line, sample, line * sample])
    fitted = smooth @ np.linalg.lstsq(smooth, kept, rcond=None)[0]
    systematic_rms = tuple(float(value) for value in np.sqrt(np.mean(fitted ** 2, axis=0)))
    return WhatAFitLeft(trusted_rms, int(np.count_nonzero(~trusted)), systematic_rms)


def where_the_references_are_seen(gcps, ref_lons, ref_lats, heights, start_time, tle, time_offset=0.0):
    """Return where the satellite sees references standing at *heights* metres above the ellipsoid.

    A point above the ellipsoid is seen along the line of sight beyond it, displaced away
    from the satellite by its height times the tangent of the satellite's zenith angle
    there. The displacement is laid off along the ground from the true position, which is
    exact to well under a metre for the few kilometres it reaches. The satellite's direction
    is taken at each gcp's own scan time, *time_offset* seconds after *start_time*.
    """
    from pyorbital.orbital import Orbital

    gcps = np.asarray(gcps, dtype=float)
    heights = np.nan_to_num(np.asarray(heights, dtype=float), nan=0.0)
    seen_at = (np.datetime64(start_time)
               + np.round((time_offset + gcps[:, 0] / 6 + gcps[:, 1] * 25e-6) * 1e9).astype("timedelta64[ns]"))
    orbit = Orbital("avhrr", line1=tle[0], line2=tle[1])
    azimuth, elevation = orbit.get_observer_look(seen_at, np.asarray(ref_lons, dtype=float),
                                                 np.asarray(ref_lats, dtype=float), heights / 1e3)
    away = heights * np.tan(np.deg2rad(90.0 - elevation))
    lons, lats, _ = geod.fwd(ref_lons, ref_lats, np.asarray(azimuth) + 180.0, away)
    return lons, lats
