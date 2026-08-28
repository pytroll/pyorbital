"""Module to compute geolocalization of a satellite scene."""

import logging
import math
import warnings
from functools import cache

import numpy as np
from pyproj import Transformer

from pyorbital import astronomy
from pyorbital.config import config
from pyorbital.orbital import Orbital

logger = logging.getLogger(__name__)

try:
    import numba as nb
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False
    logger.warning("Numba is not available, falling back to the pyproj and numpy paths.")

A = 6378.137  # WGS84 and GRS80 Equatorial radius (km)
B = 6356.752314245  # km, WGS84


@cache
def _get_transformer():
    """Get a geocentric to lon/lat transformer.

    The ellipsoid is taken from the module-level *A* and *B* constants, which the
    numba kernels read too, so every path in this module shares one declaration.
    The result is cached to avoid re-creating the PROJ context on every call.
    """
    return Transformer.from_crs(
        dict(proj="geocent", a=A * 1000, b=B * 1000),
        dict(proj="latlong", a=A * 1000, b=B * 1000))


def _ellipsoid_radius_scaling():
    """Return the scaling that maps the declared ellipsoid onto a unit sphere."""
    return np.array([[1 / A, 1 / A, 1 / B]]).T


OMEGA_EARTH = 7.2921159e-5  # Earth's rotation rate (rad/s)
_MAX_SCAN_ROWS_PER_CHUNK = 16


def geodetic_lat(point, a=A, b=B):
    """Get the Geodetic latitude of a point."""
    x, y, z = point
    r = np.sqrt(x * x + y * y)
    geoc_lat = np.arctan2(z, r)

    geod_lat = geoc_lat
    e2 = (a * a - b * b) / (a * a)
    while True:
        phi = geod_lat
        C = 1 / np.sqrt(1 - e2 * np.sin(phi)**2)
        geod_lat = np.arctan2(z + a * C * e2 * np.sin(phi), r)
        if np.allclose(geod_lat, phi):
            return geod_lat


def subpoint(query_point, a=A, b=B):
    """Get the point on the ellipsoid under the *query_point*."""
    x, y, z = query_point

    lat = geodetic_lat(query_point)
    lon = np.arctan2(y, x)
    e2_ = (a * a - b * b) / (a * a)
    n__ = a / np.sqrt(1 - e2_ * np.sin(lat)**2)
    nx_ = n__ * np.cos(lat) * np.cos(lon)
    ny_ = n__ * np.cos(lat) * np.sin(lon)
    nz_ = (1 - e2_) * n__ * np.sin(lat)

    return np.stack([nx_, ny_, nz_], axis=0)


def compute_yaw_steering(pos, vel):
    """Compute the yaw steering angle to compensate for Earth's rotation.

    Args:
        pos: Satellite position as column vector(s) in km (ECI frame).
        vel: Satellite velocity as column vector(s) in km/s (ECI frame).

    Returns:
        Yaw steering angle in radians.
    """
    r = vnorm(pos)
    v = vnorm(vel)
    lat = np.arcsin(pos[2] / r)
    return np.arctan2(OMEGA_EARTH * A * np.cos(lat), v)


def _warn_legacy_convention(what, replacement):
    """Warn that a legacy geometry convention is in use, and how to opt out."""
    warnings.warn(
        f"pyorbital is using the legacy {what}. It is the default so that products "
        "generated with earlier versions remain reproducible, but validation against "
        f"reference geolocation shows it is less accurate; pass {replacement} for the "
        "corrected computation. The default will change in a future release.",
        DeprecationWarning,
        stacklevel=3,
    )


def _resolve_convention(value, config_key, description, replacement):
    """Resolve a geometry convention from the argument, then the config, then legacy.

    Downstream libraries call into the geolocation without forwarding these
    arguments, so the process-wide config is what lets an archive reprocessing
    select a convention.  Both an explicit argument and a config entry count as a
    deliberate choice and stay silent; only falling through to the legacy default
    warns.
    """
    if value is not None:
        return value
    configured = config.get(config_key, None)
    if configured is not None:
        return configured
    _warn_legacy_convention(description, replacement)
    return "legacy"


def _local_frame(pos, vel, nadir_convention=None):
    """Compute the satellite's local orbital reference frame.

    Returns (nadir, along_track, cross_track) as unit column vectors.

    *nadir_convention* selects the "down" direction:

    ``"legacy"``
        The historical direction, the normalised position of the ellipsoid
        subpoint of the antipode.  Used when *nadir_convention* is left
        unspecified, so that products already generated with pyorbital reproduce
        unchanged.  Relying on that default emits a :exc:`DeprecationWarning`,
        because the convention is measurably less accurate; requesting
        ``"legacy"`` explicitly is treated as an informed choice and is silent.
    ``"geocentric"``
        Straight at the centre of the Earth.  Validation of FY-3 MERSI against
        reference geolocation selects this one.
    ``"geodetic"``
        Along the ellipsoid normal, from the satellite to ``subpoint(pos)``.

    The nadir is re-orthogonalised against ``along_track`` via Gram-Schmidt to
    maintain the orthonormality required by the broadcast rotation in
    :class:`ScanGeometry`.
    """
    nadir_convention = _resolve_convention(
        nadir_convention, "nadir_convention", "nadir convention", "nadir_convention='geocentric'"
    )
    along_track = vel / vnorm(vel)
    if nadir_convention == "legacy":
        # Reproduce released pyorbital exactly, which means *not* re-orthogonalising:
        # the Gram-Schmidt below is itself a behavioural change, worth up to ~2.7 km
        # once a pitch bias is involved.  Callers asking for the legacy convention
        # are asking for the released numbers, so they get the released frame.
        nadir = subpoint(-pos)
        nadir = nadir / vnorm(nadir)
    else:
        if nadir_convention == "geocentric":
            nadir = -pos / vnorm(pos)
        else:
            nadir = subpoint(pos) - pos
            nadir = nadir / vnorm(nadir)
        # The rotations need an orthonormal triad, but the requested nadir is the
        # whole point of the convention, so the along-track axis absorbs the
        # correction rather than the nadir.  A real orbit always has some radial
        # velocity -- the J2 short-period oscillation alone gives a flight-path
        # angle of order 1e-3 rad -- and tilting the nadir by that would put ~1 km
        # on the ground, flipping sign between the ascending and descending legs.
        cross_track = np.cross(nadir, vel, 0, 0, 0)
        cross_track = cross_track / vnorm(cross_track)
        along_track = np.cross(cross_track, nadir, 0, 0, 0)
        along_track = along_track / vnorm(along_track)
        return nadir, along_track, cross_track
    cross_track = np.cross(nadir, vel, 0, 0, 0)
    cross_track = cross_track / vnorm(cross_track)
    return nadir, along_track, cross_track


def _effective_yaw(yaw, yaw_steering, pos, vel, fovs_shape):
    """Compute the effective yaw angle, optionally adding the yaw-steering term."""
    if yaw_steering:
        yaw = yaw + compute_yaw_steering(pos, vel)
    if np.shape(yaw):
        yaw_arr = np.asarray(yaw)
        if yaw_arr.ndim == 1 and len(fovs_shape) == 2:
            yaw_arr = yaw_arr[:, np.newaxis]   # (M,) → (M, 1) broadcasts to (M, N)
        yaw = np.broadcast_to(yaw_arr, fovs_shape)
    return yaw


def _cross3(a, b):
    """Cross product of two 3-vectors stored along axis 0.

    Works for any pair of broadcastable shapes whose first dimension is 3.
    """
    out = np.empty(np.broadcast_shapes(a.shape, b.shape), dtype=np.result_type(a, b))
    out[0] = a[1] * b[2] - a[2] * b[1]
    out[1] = a[2] * b[0] - a[0] * b[2]
    out[2] = a[0] * b[1] - a[1] * b[0]
    return out


def _rodrigues(vector, axis, angle):
    """Rotate *vector* around *axis* by *angle* using the Rodrigues formula.

    **Sign convention** — identical to :func:`qrotate`: the rotation is
    *clockwise* when viewed from the positive end of *axis*.  Concretely,
    rotating the x-axis around the z-axis by +90° yields ``[0, -1, 0]``,
    not ``[0, +1, 0]``.  This matches the (intentional) sign of the existing
    quaternion path; instrument scan-angle arrays are defined with the
    appropriate sign to produce the correct swath geometry.

    Supports broadcasting: *vector* and *axis* may have fewer trailing
    dimensions than *angle*.  The two handled cases are:

    * ``vector`` / ``axis`` shape ``(3, M)`` and ``angle`` shape ``(M, N)``:
      returns ``(3, M, N)`` — per-scan axis applied to per-pixel angles.
    * ``vector`` shape ``(3, M, N)``, ``axis`` ``(3, M)``, ``angle``
      ``(M, N)``: returns ``(3, M, N)`` — used for the second rotation pass.

    For flat (2-D) inputs the behaviour is identical to :func:`qrotate`.
    """
    n = axis / vnorm(axis)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    if vector.ndim == 2 and np.ndim(angle) == 2:
        # Per-scan vector/axis, per-pixel angle → broadcast to (3, M, N)
        n3 = n[:, :, np.newaxis]
        v3 = vector[:, :, np.newaxis]
        c3 = cos_a[np.newaxis]
        s3 = sin_a[np.newaxis]
        ndotv = (n * vector).sum(0)
        ncrossv = _cross3(n, vector)
        return (v3 * c3
                - ncrossv[:, :, np.newaxis] * s3
                + n3 * ndotv[np.newaxis, :, np.newaxis] * (1 - c3))

    if vector.ndim == 3 and axis.ndim == 2 and np.ndim(angle) == 2:
        # vector (3, M, N); per-scan axis (3, M); per-pixel angle (M, N)
        n3 = n[:, :, np.newaxis]
        c3 = cos_a[np.newaxis]
        s3 = sin_a[np.newaxis]
        ndotv = (n3 * vector).sum(0, keepdims=True)
        ncrossv = _cross3(n3, vector)
        n3_onec3 = n3 * (1.0 - c3)
        return vector * c3 - ncrossv * s3 + n3_onec3 * ndotv

    # Flat path: all inputs share the same trailing dimensions
    ndotv = (n * vector).sum(0)
    ncrossv = _cross3(n, vector)
    return vector * cos_a - ncrossv * sin_a + n * ndotv * (1 - cos_a)


class ScanGeometry(object):
    """Description of the geometry of an instrument.

    *fovs* is the x and y viewing angles of the instrument. y is zero if the we
    talk about scanlines of course. *times* is the time of viewing of each
    angle relative to the start of the scanning, so it should have the same
    size as the *fovs*. *attitude* is the attitude correction to apply.
    """

    def __init__(self, fovs, times, attitude=(0, 0, 0), lines_per_scan=1):
        """Initialize the class."""
        self.fovs = np.array(fovs)
        self.lines_per_scan = int(lines_per_scan)
        try:
            # assuming seconds
            self._times = np.asanyarray(times) * np.timedelta64(1000000000, "ns")
        except TypeError:
            self._times = np.asanyarray(times).astype("timedelta64[ns]")
        self.attitude = attitude

    def vectors(self, pos, vel, roll=0.0, pitch=0.0, yaw=0.0, yaw_steering=False,
                nadir_convention=None, rotation_order=None):
        """Get unit vectors pointing to the different pixels.

        *pos* and *vel* are column vectors, or matrices of column
        vectors. Returns vectors as stacked rows.

        If *yaw_steering* is True, the yaw angle is computed from the
        satellite position and velocity to compensate for Earth's rotation.
        This is added to any explicit *yaw* value.

        When fovs is 3-D (shape ``(2, N_scans, N_pixels)``) and *pos*/*vel*
        are per-scan (shape ``(3, N_scans)``), a broadcasting Rodrigues
        rotation is used to avoid repeating the local frame computation for
        every pixel.  The result is then flattened to ``(3, N_total)``.
        """
        nadir_convention = _resolve_convention(
            nadir_convention, "nadir_convention", "nadir convention",
            "nadir_convention='geocentric'"
        )
        per_scan_state = (self.fovs.ndim == 3 and pos.ndim == 2
                          and pos.shape[1] == self.fovs.shape[1])
        if per_scan_state and nadir_convention != "legacy":
            return self._vectors_broadcast(pos, vel, roll, pitch, yaw, yaw_steering,
                                           nadir_convention, rotation_order)
        if per_scan_state:
            # The broadcast rotation simplifies using nadir ⊥ along_track, which the
            # legacy frame deliberately does not satisfy, so legacy takes the general
            # path instead -- which needs the orbit state repeated for every pixel.
            pixels_per_row = self.fovs.shape[2]
            pos = np.repeat(pos, pixels_per_row, axis=1)
            vel = np.repeat(vel, pixels_per_row, axis=1)
        fovs = self.fovs.reshape(2, -1)
        nadir, along_track, cross_track = _local_frame(pos, vel, nadir_convention)
        effective_yaw = _effective_yaw(yaw, yaw_steering, pos, vel, self.fovs[0].shape)
        if not np.any(fovs[1] + pitch):
            # with no along-track angle the two orders are identical, so the
            # choice is irrelevant and not worth deprecating at the caller
            rotated = qrotate(nadir, along_track, fovs[0] + roll)
            return qrotate(rotated, nadir, effective_yaw)
        rotation_order = _resolve_convention(
            rotation_order, "rotation_order", "rotation order", "rotation_order='pitch_first'"
        )
        if rotation_order == "legacy":
            # as released: roll (about along-track) first, then pitch
            rotated = qrotate(nadir, along_track, fovs[0] + roll)
            rotated = qrotate(rotated, cross_track, fovs[1] + pitch)
        else:
            rotated = qrotate(nadir, cross_track, fovs[1] + pitch)
            rotated = qrotate(rotated, along_track, fovs[0] + roll)
        return qrotate(rotated, nadir, effective_yaw)

    def _vectors_broadcast(self, pos, vel, roll, pitch, yaw, yaw_steering,
                           nadir_convention=None, rotation_order=None):
        """Broadcast rotation for per-scan pos/vel and per-pixel fovs.

        Two structural invariants of all 3-D fov scan geometries are exploited:

        1. **Cross-track angles** (``fovs[0]``) are the same for every scan line
           — they only depend on the pixel index.  Trig is therefore computed
           on ``N_pixels`` values instead of ``N_scans * N_pixels``.

        2. **Along-track angles** (``fovs[1]``) are constant across all pixels
           within a scan line — they only depend on the detector-row index.
           Trig is computed on ``N_scans`` values instead of
           ``N_scans * N_pixels``.

        3. **nadir ⊥ cross_track** (guaranteed by :func:`_local_frame`).
           The along-track Rodrigues rotation therefore simplifies to two
           multiply-add operations:
           ``r1 = nadir * cos(α) - along_track * sin(α)``
        """
        nadir, along_track, cross_track = _local_frame(pos, vel, nadir_convention)
        effective_yaw = _effective_yaw(yaw, yaw_steering, pos, vel, self.fovs[0].shape)

        nadir_3d = nadir[:, :, np.newaxis]                      # (3, M, 1)
        at_3d = along_track[:, :, np.newaxis]                   # (3, M, 1)
        ct_3d = cross_track[:, :, np.newaxis]                   # (3, M, 1)

        # Cross-track angles are identical for all scan lines: compute trig once.
        cross_angles = (self.fovs[0] + roll)[0:1, :]           # (1, N)
        cos_c = np.cos(cross_angles)                            # (1, N)
        sin_c = np.sin(cross_angles)                            # (1, N)

        # Along-track angles are constant across pixels: compute trig once per row.
        along_angles = (self.fovs[1] + pitch)[:, 0:1]          # (M, 1)

        if np.any(along_angles):
            rotation_order = _resolve_convention(
                rotation_order, "rotation_order", "rotation order",
                "rotation_order='pitch_first'"
            )
            cos_a = np.cos(along_angles)[np.newaxis]            # (1, M, 1)
            sin_a = np.sin(along_angles)[np.newaxis]            # (1, M, 1)
            if rotation_order == "legacy":
                # Cross-track first (nadir ⊥ along_track makes that the simple one),
                # then along-track.  Closed form, matching the flat path and kernel:
                # r = nadir*cos(θ)*cos(α) - along_track*cos(θ)*sin(α)
                #     + cross_track*sin(θ)
                rotated = (nadir_3d * cos_c[np.newaxis] * cos_a
                           - at_3d * cos_c[np.newaxis] * sin_a
                           + ct_3d * sin_c[np.newaxis])
            else:
                # --- Step 1: along-track (simplified Rodrigues around cross_track) ---
                # Exploit nadir ⊥ cross_track: r1 = nadir*cos(α) - along_track*sin(α)
                r1 = nadir_3d * cos_a - at_3d * sin_a              # (3, M, 1)

                # --- Step 2: cross-track (compact full-Rodrigues around along_track) ---
                # r2 = r1*cos(θ) + cross_track*cos(α)*sin(θ)
                #      - along_track*sin(α)*(1-cos(θ))
                one_m_cos_c = 1.0 - cos_c                          # (1, N)
                rotated = (r1 * cos_c[np.newaxis]
                           + ct_3d * cos_a * sin_c[np.newaxis]
                           - at_3d * sin_a * one_m_cos_c[np.newaxis])
        else:
            # No along-track offset: simplified cross-track only (nadir ⊥ along_track).
            rotated = nadir_3d * cos_c[np.newaxis] + ct_3d * sin_c[np.newaxis]

        # --- Yaw rotation (only if non-zero) ---
        if np.shape(effective_yaw):
            effective_yaw = effective_yaw.reshape(self.fovs[0].shape)
        if np.any(effective_yaw):
            rotated = _rodrigues(rotated, nadir, np.broadcast_to(effective_yaw, self.fovs[0].shape))

        return rotated.reshape(3, -1)                           # (3, N_total)

    def times(self, start_of_scan):
        """Return an array with the times of each scan line."""
        try:
            return np.array(self._times) + np.datetime64(start_of_scan)
        except ValueError:
            return np.array(self._times) + start_of_scan


class Quaternion(object):
    """Some class, that I don't know what is doing..."""

    def __init__(self, scalar, vector):
        """Initialize the class."""
        self.__x, self.__y, self.__z = vector.reshape((3, -1))
        self.__w = scalar.ravel()

    def rotation_matrix(self):
        """Get the rotation matrix."""
        x, y, z, w = self.__x, self.__y, self.__z, self.__w
        zero = np.zeros_like(x)
        return np.array(
            ((w**2 + x**2 - y**2 - z**2,
              2 * x * y + 2 * z * w,
              2 * x * z - 2 * y * w,
              zero),
             (2 * x * y - 2 * z * w,
              w**2 - x**2 + y**2 - z**2,
              2 * y * z + 2 * x * w,
              zero),
             (2 * x * z + 2 * y * w,
              2 * y * z - 2 * x * w,
              w**2 - x**2 - y**2 + z**2,
              zero),
             (zero, zero, zero, w**2 + x**2 + y**2 + z**2)))


def qrotate(vector, axis, angle):
    """Rotate *vector* around *axis* by *angle* (in radians).

    *vector* is a matrix of column vectors, as is *axis*.
    This function uses quaternion rotation.

    **Sign convention**: the rotation is *clockwise* when viewed from the
    positive end of *axis* (i.e. the right-hand-rule angle is ``-angle``).
    This is equivalent to applying the conjugate quaternion q†vq instead of
    the standard qvq†. All scan-angle arrays in the instrument definitions
    are negated accordingly so that positive scan angles map to the
    right-hand side of the swath.

    Use :func:`_rodrigues` for new code, as it encodes the same convention but
    supports broadcasting across scan lines.
    """
    n_axis = axis / vnorm(axis)
    sin_angle = np.expand_dims(np.sin(angle / 2), 0)
    if np.ndim(n_axis) == 1:
        n_axis = np.expand_dims(n_axis, 1)
        p__ = np.dot(n_axis, sin_angle)[:, np.newaxis]
    else:
        p__ = n_axis * sin_angle

    q__ = Quaternion(np.cos(angle / 2), p__)
    shape = vector.shape
    return np.einsum("kj, ikj->ij",
                     vector.reshape((3, -1)),
                     q__.rotation_matrix()[:3, :3]).reshape(shape)


if _HAS_NUMBA:
    @nb.njit(parallel=True, cache=True)
    def _numba_ecef_to_lonlatalt(x, y, z):
        """Convert flat ECEF arrays (km) to lon (deg), lat (deg), alt (m).

        Implements Bowring's one-iteration geodetic formula in a fully parallel
        numba JIT kernel.  Accuracy: < 1e-10 deg for lon/lat, < 1e-4 m for
        altitude — identical to pyproj's PROJ implementation.

        Args:
            x: 1-D float64 array of ECEF x coordinates in kilometres.
            y: 1-D float64 array of ECEF y coordinates in kilometres.
            z: 1-D float64 array of ECEF z coordinates in kilometres.

        Returns:
            Tuple ``(lon_deg, lat_deg, alt_m)`` as 1-D float64 arrays.
        """
        a = A
        b = B
        _e2 = 1.0 - (b / a) ** 2
        _ep2 = (a / b) ** 2 - 1.0
        _r2d = 180.0 / math.pi
        _sin_89_9 = 0.9998476951563913  # sin(89.9°), pole threshold

        n = len(x)
        lon_out = np.empty(n)
        lat_out = np.empty(n)
        alt_out = np.empty(n)

        for i in nb.prange(n):
            xi, yi, zi = x[i], y[i], z[i]
            lon_out[i] = math.atan2(yi, xi) * _r2d

            p = math.sqrt(xi * xi + yi * yi)
            za = zi * a
            pb = p * b
            r_t = math.sqrt(za * za + pb * pb)
            sin_t = za / r_t
            cos_t = pb / r_t

            num = zi + _ep2 * b * sin_t * sin_t * sin_t
            den = p - _e2 * a * cos_t * cos_t * cos_t
            lat = math.atan2(num, den)
            lat_out[i] = lat * _r2d

            sin_lat = math.sin(lat)
            cos_lat = math.cos(lat)
            N = a / math.sqrt(1.0 - _e2 * sin_lat * sin_lat)
            if abs(sin_lat) > _sin_89_9:
                alt_out[i] = (abs(zi) / abs(sin_lat) - N * (1.0 - _e2)) * 1000.0
            else:
                alt_out[i] = (p / cos_lat - N) * 1000.0

        return lon_out, lat_out, alt_out

    @nb.njit(parallel=True, cache=True)
    def _fused_geolocate_numba(pos, nadir, cross_track, cos_c, sin_c, along_angles,
                               pitch_first,
                               yaw_angles, gmst_scan, lon_out, lat_out, alt_out):
        """Fused rotation + ellipsoid intersection + geodetic in one parallel kernel.

        Outer ``prange`` over M scan lines (one per CPU thread); inner sequential
        loop over N pixels.  Eliminates the large ``(3, M, N)`` intermediate
        arrays produced by the numpy pipeline, reducing memory pressure by ~3×.

        Sign convention: clockwise Rodrigues (negative sin term) matching
        :func:`_rodrigues` and :func:`qrotate`.

        Three rotations are applied in order:
        1. Along-track (simplified Rodrigues, exploits nadir ⊥ cross_track).
        2. Cross-track (compact full-Rodrigues around along_track = cross_track × nadir).
        3. Yaw Rodrigues around nadir (skipped when ``yaw_angles[m] ≈ 0``).

        Args:
            pos: Satellite ECEF positions (km), shape ``(3, M)``.
            nadir: Pre-computed geodetic nadir unit vectors from
                :func:`_local_frame`, shape ``(3, M)``.
            cross_track: Pre-computed cross-track unit vectors from
                :func:`_local_frame`, shape ``(3, M)``.
            pitch_first: Whether to apply pitch before roll; ``False`` selects
                the released roll-then-pitch order.
            cos_c: cos of (cross-track scan angle + roll), shape ``(N,)``.
            sin_c: sin of (cross-track scan angle + roll), shape ``(N,)``.
            along_angles: Along-track detector angle + pitch per scan line,
                shape ``(M,)`` in radians.
            yaw_angles: Effective yaw per scan line (explicit yaw + optional
                steering), shape ``(M,)`` in radians.
            gmst_scan: Greenwich Mean Sidereal Time per scan line, shape
                ``(M,)`` in radians.
            lon_out: Pre-allocated output array of longitudes in degrees,
                shape ``(M * N,)``.
            lat_out: Pre-allocated output array of latitudes in degrees,
                shape ``(M * N,)``.
            alt_out: Pre-allocated output array of altitudes in metres,
                shape ``(M * N,)``.
        """
        _a = A
        _b = B
        _e2 = 1.0 - (_b / _a) ** 2
        _ep2 = (_a / _b) ** 2 - 1.0
        _ra = 1.0 / _a
        _rb = 1.0 / _b
        _r2d = 180.0 / math.pi
        _pi = math.pi
        _twopi = 2.0 * _pi
        _sin_89_9 = 0.9998476951563913

        M = pos.shape[1]
        N = len(cos_c)

        for m in nb.prange(M):
            px, py, pz = pos[0, m], pos[1, m], pos[2, m]

            nx, ny, nz = nadir[0, m], nadir[1, m], nadir[2, m]
            ctx, cty, ctz = cross_track[0, m], cross_track[1, m], cross_track[2, m]

            # Along-track trig — computed once per scan line
            alpha = along_angles[m]
            cos_a = math.cos(alpha)
            sin_a = math.sin(alpha)

            # along_track = cross_track × nadir  (pre-computed per row)
            atx = cty * nz - ctz * ny
            aty = ctz * nx - ctx * nz
            atz = ctx * ny - cty * nx

            # --- Step 1: along-track (simplified Rodrigues, nadir ⊥ cross_track) ---
            # r1 = nadir*cos(α) - along_track*sin(α)  (pre-computed per row)
            r1x = nx * cos_a - atx * sin_a
            r1y = ny * cos_a - aty * sin_a
            r1z = nz * cos_a - atz * sin_a

            gmst_m = gmst_scan[m]

            for k in range(N):
                cc, sc = cos_c[k], sin_c[k]
                if pitch_first:
                    # --- Step 2: cross-track (compact full-Rodrigues around along_track) ---
                    # r2 = r1*cos(θ) + cross_track*cos(α)*sin(θ)
                    #      - along_track*sin(α)*(1-cos(θ))
                    one_m_cc = 1.0 - cc
                    rx = r1x * cc + ctx * cos_a * sc - atx * sin_a * one_m_cc
                    ry = r1y * cc + cty * cos_a * sc - aty * sin_a * one_m_cc
                    rz = r1z * cc + ctz * cos_a * sc - atz * sin_a * one_m_cc
                else:
                    # Legacy order: cross-track first (nadir ⊥ along_track, so that
                    # rotation is the simple one), then along-track.  Closed form:
                    # r = nadir*cos(θ)*cos(α) - along_track*cos(θ)*sin(α)
                    #     + cross_track*sin(θ)
                    rx = nx * cc * cos_a - atx * cc * sin_a + ctx * sc
                    ry = ny * cc * cos_a - aty * cc * sin_a + cty * sc
                    rz = nz * cc * cos_a - atz * cc * sin_a + ctz * sc

                # --- 3rd rotation: CW Rodrigues around nadir by yaw angle ---
                yaw_m = yaw_angles[m]
                cos_y = math.cos(yaw_m)
                sin_y = math.sin(yaw_m)
                if abs(sin_y) > 1e-15:
                    ndotr = nx * rx + ny * ry + nz * rz
                    nxrx = ny * rz - nz * ry
                    nxry = nz * rx - nx * rz
                    nxrz = nx * ry - ny * rx
                    one_m_cos_y = 1.0 - cos_y
                    rx = rx * cos_y - nxrx * sin_y + nx * ndotr * one_m_cos_y
                    ry = ry * cos_y - nxry * sin_y + ny * ndotr * one_m_cos_y
                    rz = rz * cos_y - nxrz * sin_y + nz * ndotr * one_m_cos_y

                # --- Ellipsoid intersection (WGS-84) ---
                xr = rx * _ra
                yr = ry * _ra
                zr = rz * _rb
                cxr = -px * _ra
                cyr = -py * _ra
                czr = -pz * _rb
                ldotc = xr * cxr + yr * cyr + zr * czr
                lsq = xr * xr + yr * yr + zr * zr
                csq = cxr * cxr + cyr * cyr + czr * czr
                disc = ldotc * ldotc - csq * lsq + lsq
                if disc < 0.0:
                    disc = 0.0
                d1 = (ldotc - math.sqrt(disc)) / lsq
                ex = rx * d1 + px
                ey = ry * d1 + py
                ez = rz * d1 + pz

                # --- Bowring geodetic: ECEF (km) → lon/lat/alt ---
                lon = math.atan2(ey, ex) - gmst_m
                if lon < -_pi:
                    lon += _twopi
                elif lon > _pi:
                    lon -= _twopi

                p = math.sqrt(ex * ex + ey * ey)
                za = ez * _a
                pb = p * _b
                r_t = math.sqrt(za * za + pb * pb)
                sin_t = za / r_t
                cos_t = pb / r_t
                num = ez + _ep2 * _b * sin_t * sin_t * sin_t
                den = p - _e2 * _a * cos_t * cos_t * cos_t
                lat = math.atan2(num, den)

                sin_lat = math.sin(lat)
                cos_lat = math.cos(lat)
                N_val = _a / math.sqrt(1.0 - _e2 * sin_lat * sin_lat)
                if abs(sin_lat) > _sin_89_9:
                    alt = (abs(ez) / abs(sin_lat) - N_val * (1.0 - _e2)) * 1000.0
                else:
                    alt = (p / cos_lat - N_val) * 1000.0

                idx = m * N + k
                lon_out[idx] = lon * _r2d
                lat_out[idx] = lat * _r2d
                alt_out[idx] = alt

else:
    def _numba_ecef_to_lonlatalt(x, y, z):  # type: ignore[misc]
        """Stub: numba not available — callers should use pyproj instead."""
        raise RuntimeError("numba is not installed")

    def _fused_geolocate_numba(pos, nadir, cross_track, cos_c, sin_c, along_angles,  # type: ignore[misc]
                               pitch_first, yaw_angles, gmst_scan, lon_out, lat_out, alt_out):
        """Stub: numba not available."""
        raise RuntimeError("numba is not installed")


def get_lonlatalt(pos, utc_time):
    """Calculate sublon, sublat and altitude of satellite, considering the earth an ellipsoid."""
    if _HAS_NUMBA:
        lon, lat, alt = _numba_ecef_to_lonlatalt(pos[0], pos[1], pos[2])
    else:
        lon, lat, alt = _get_transformer().transform(*(pos * 1000))
    lon = lon - np.rad2deg(_gmst_per_pixel(utc_time))
    lon = np.where(lon < -180, lon + 360, lon)
    return lon, lat, alt


def get_sensor_angles(orb, utc_time, lon, lat, alt=0.0):
    """Return sensor zenith and azimuth angles for ground coordinates.

    Coordinates are degrees and ground altitude is kilometres above the geoid.
    The orbit provider must implement ``get_position(utc_time, normalize=False)``.
    """
    from pyorbital.orbital import get_observer_look

    if isinstance(orb, (list, tuple)):
        orb = Orbital("mysatellite", line1=orb[0], line2=orb[1])
    pos, _ = orb.get_position(utc_time, normalize=False)
    pos = np.asarray(pos).reshape(3, -1)
    sat_lon, sat_lat, sat_alt = get_lonlatalt(pos, utc_time)
    azimuth, elevation = get_observer_look(sat_lon, sat_lat, sat_alt / 1000.0, utc_time, lon, lat, alt)
    return 90.0 - elevation, azimuth % 360.0


def geolocate(orb, sgeom, times, rpy=(0.0, 0.0, 0.0), yaw_steering=False,
              nadir_convention=None, rotation_order=None):
    """Compute (lon, lat, alt) for every pixel in one shot.

    When numba is available **and** the scan geometry has 3-D fovs (any
    :class:`~pyorbital.geoloc_instrument_definitions.MultiLineSweepbroomScan`
    or equivalent), a single fused parallel JIT kernel handles rotation
    (including optional yaw steering) + ellipsoid intersection + geodetic
    conversion with no large intermediate arrays.  Otherwise falls back to
    :func:`compute_pixels` + :func:`get_lonlatalt`.

    Args:
        orb: :class:`~pyorbital.orbital.Orbital` instance *or* ``(tle1, tle2)``
            tuple.
        sgeom: :class:`~pyorbital.geoloc.ScanGeometry` instance.
        times: Per-pixel UTC times — 2-D array ``(N_rows, N_pixels)`` returned
            by :meth:`~pyorbital.geoloc.ScanGeometry.times`.
        rpy: ``(roll, pitch, yaw)`` attitude corrections in radians, or one
            such row per instrument scan.
        yaw_steering: If ``True``, compute yaw from orbit geometry to
            counteract Earth rotation.
        nadir_convention: Which "down" direction the local frame uses; see
            :func:`_local_frame`.  Left unspecified it keeps the legacy
            convention and emits a :exc:`DeprecationWarning`.
        rotation_order: ``"legacy"`` applies roll before pitch, as released;
            ``"pitch_first"`` applies pitch before roll so the along-track
            component does not depend on the scan angle.  Left unspecified it
            keeps the legacy order and emits a :exc:`DeprecationWarning`.

    Returns:
        Tuple ``(lon_deg, lat_deg, alt_m)`` as flat 1-D arrays.
    """
    rpy_array = np.asarray(rpy)
    if rpy_array.ndim == 2:
        return _geolocate_per_scan_attitude(orb, sgeom, times, rpy_array, yaw_steering,
                                            nadir_convention, rotation_order)
    roll, pitch, yaw = rpy

    times_arr = np.asanyarray(times)
    if _requires_per_pixel_orbit(sgeom, times_arr):
        return _geolocate_with_per_pixel_orbit(orb, sgeom, times_arr, rpy, yaw_steering,
                                               nadir_convention, rotation_order)

    if _HAS_NUMBA and sgeom.fovs.ndim == 3:
        return _geolocate_fused(orb, sgeom, times, roll, pitch, yaw, yaw_steering,
                                nadir_convention, rotation_order)

    pixels = compute_pixels(orb, sgeom, times, rpy, yaw_steering, nadir_convention,
                            rotation_order)
    return get_lonlatalt(pixels, times)


def _geolocate_per_scan_attitude(orb, sgeom, times, attitudes, yaw_steering,
                                 nadir_convention=None, rotation_order=None):
    """Geolocate compact scans carrying independent attitude corrections."""
    rows_per_scan = sgeom.lines_per_scan
    parts = []
    for scan_index, attitude in enumerate(attitudes):
        row_start = scan_index * rows_per_scan
        row_stop = row_start + rows_per_scan
        scan_times = np.asanyarray(times)[row_start:row_stop]
        scan_geometry = ScanGeometry(
            sgeom.fovs[:, row_start:row_stop],
            scan_times - scan_times[0, 0],
            lines_per_scan=rows_per_scan,
        )
        parts.append(geolocate(orb, scan_geometry, scan_times, tuple(attitude), yaw_steering,
                               nadir_convention, rotation_order))
    return tuple(np.concatenate(values) for values in zip(*parts))


def _requires_per_pixel_orbit(sgeom, times):
    """Return whether pixels within a compact scan have distinct acquisition times."""
    return sgeom.fovs.ndim == 3 and times.ndim == 2 and np.any(times != times[:, :1])


def _geolocate_with_per_pixel_orbit(orb, sgeom, times, rpy, yaw_steering,
                                    nadir_convention=None, rotation_order=None):
    """Geolocate a compact scan without discarding its intra-scan orbital motion."""
    results = []
    for start in range(0, times.shape[0], _MAX_SCAN_ROWS_PER_CHUNK):
        stop = min(start + _MAX_SCAN_ROWS_PER_CHUNK, times.shape[0])
        chunk_geometry = ScanGeometry(sgeom.fovs[:, start:stop], times[start:stop] - times[start, 0])
        results.append(_geolocate_scan_chunk(orb, chunk_geometry, times[start:stop], rpy,
                                             yaw_steering, nadir_convention, rotation_order))
    return tuple(np.concatenate(parts) for parts in zip(*results))


def _geolocate_scan_chunk(orb, sgeom, times, rpy, yaw_steering, nadir_convention=None,
                          rotation_order=None):
    """Geolocate a bounded group of compact scan rows."""
    pixel_times = times.reshape(-1)
    flat_fovs = sgeom.fovs.reshape(2, -1)
    pos, vel = _interpolate_scan_endpoint_states(orb, times)

    # The fused kernel covers the configuration FY-3 MERSI uses. Anything else
    # — yaw steering, the legacy conventions — falls through to the array path,
    # which remains the reference implementation.
    if (_HAS_NUMBA and not yaw_steering
            and nadir_convention == "geocentric" and rotation_order == "pitch_first"):
        roll, pitch, yaw = rpy
        return _geocentric_scan_chunk_numba(
            np.ascontiguousarray(flat_fovs),
            np.ascontiguousarray(pos),
            np.ascontiguousarray(vel),
            float(roll), float(pitch), float(yaw),
            np.ascontiguousarray(_gmst_per_pixel(pixel_times)),
            bool(np.any(flat_fovs[1] + pitch)),
        )

    flat_geometry = ScanGeometry(flat_fovs, pixel_times - pixel_times[0])
    vectors = flat_geometry.vectors(pos, vel, *rpy, yaw_steering=yaw_steering,
                                    nadir_convention=nadir_convention,
                                    rotation_order=rotation_order)
    radius = _ellipsoid_radius_scaling()
    pixels = _ellipsoid_intersection(vectors, pos, radius)
    return get_lonlatalt(pixels, pixel_times)


if _HAS_NUMBA:
    @nb.njit(inline="always", cache=True)
    def _qrot_scalar(vx, vy, vz, ax, ay, az, angle):
        """Rotate one vector about one axis, reproducing :func:`qrotate`.

        Mirrors the quaternion -> rotation-matrix -> contraction sequence,
        including the axis re-normalisation qrotate always performs, so the
        kernel agrees with the array path to the last bits rather than merely
        to an algebraic identity.
        """
        n = math.sqrt(ax * ax + ay * ay + az * az)
        nx, ny, nz = ax / n, ay / n, az / n
        sin_half = math.sin(angle / 2)
        w = math.cos(angle / 2)
        x, y, z = nx * sin_half, ny * sin_half, nz * sin_half
        r00 = w * w + x * x - y * y - z * z
        r01 = 2 * x * y + 2 * z * w
        r02 = 2 * x * z - 2 * y * w
        r10 = 2 * x * y - 2 * z * w
        r11 = w * w - x * x + y * y - z * z
        r12 = 2 * y * z + 2 * x * w
        r20 = 2 * x * z + 2 * y * w
        r21 = 2 * y * z - 2 * x * w
        r22 = w * w - x * x - y * y + z * z
        return (vx * r00 + vy * r01 + vz * r02,
                vx * r10 + vy * r11 + vz * r12,
                vx * r20 + vy * r21 + vz * r22)

    @nb.njit(parallel=True, cache=True)
    def _geocentric_scan_chunk_numba(fovs, pos, vel, roll, pitch, yaw, gmst, has_along_track):
        """Fused per-pixel-orbit geolocation for the geocentric nadir convention.

        One pass over pixels doing local frame, rotations, ellipsoid
        intersection and the geodetic conversion, with no large intermediates.
        The satellite state is per-pixel, so each scan's intra-scan motion is
        preserved exactly as in the array path.
        """
        a = A
        b = B
        bow_e2 = 1.0 - (b / a) ** 2
        bow_ep2 = (a / b) ** 2 - 1.0
        r2d = 180.0 / math.pi
        inv_a, inv_b = 1.0 / a, 1.0 / b

        n = fovs.shape[1]
        lon_out = np.empty(n)
        lat_out = np.empty(n)
        alt_out = np.empty(n)

        for i in nb.prange(n):
            px, py, pz = pos[0, i], pos[1, i], pos[2, i]
            vx, vy, vz = vel[0, i], vel[1, i], vel[2, i]

            # Geocentric nadir: straight at the centre of the Earth.  It is kept
            # exact and the along-track axis absorbs the orthogonalisation, as in
            # _local_frame -- orthogonalising the *nadir* against the velocity
            # instead would tilt it by the flight-path angle (~1e-3 rad from the
            # J2 short-period oscillation alone, ~1 km on the ground, flipping
            # sign between the ascending and descending legs).
            pnorm = math.sqrt(px * px + py * py + pz * pz)
            nad_x, nad_y, nad_z = -px / pnorm, -py / pnorm, -pz / pnorm

            ct_x = nad_y * vz - nad_z * vy
            ct_y = nad_z * vx - nad_x * vz
            ct_z = nad_x * vy - nad_y * vx
            cn = math.sqrt(ct_x * ct_x + ct_y * ct_y + ct_z * ct_z)
            ct_x, ct_y, ct_z = ct_x / cn, ct_y / cn, ct_z / cn

            at_x = ct_y * nad_z - ct_z * nad_y
            at_y = ct_z * nad_x - ct_x * nad_z
            at_z = ct_x * nad_y - ct_y * nad_x
            an = math.sqrt(at_x * at_x + at_y * at_y + at_z * at_z)
            at_x, at_y, at_z = at_x / an, at_y / an, at_z / an

            # pitch_first: about cross-track, then roll about along-track.
            if has_along_track:
                rx, ry, rz = _qrot_scalar(nad_x, nad_y, nad_z, ct_x, ct_y, ct_z, fovs[1, i] + pitch)
                rx, ry, rz = _qrot_scalar(rx, ry, rz, at_x, at_y, at_z, fovs[0, i] + roll)
            else:
                rx, ry, rz = _qrot_scalar(nad_x, nad_y, nad_z, at_x, at_y, at_z, fovs[0, i] + roll)
            rx, ry, rz = _qrot_scalar(rx, ry, rz, nad_x, nad_y, nad_z, yaw)

            xr0, xr1, xr2 = rx * inv_a, ry * inv_a, rz * inv_b
            cr0, cr1, cr2 = -px * inv_a, -py * inv_a, -pz * inv_b
            ldotc = xr0 * cr0 + xr1 * cr1 + xr2 * cr2
            lsq = xr0 * xr0 + xr1 * xr1 + xr2 * xr2
            csq = cr0 * cr0 + cr1 * cr1 + cr2 * cr2
            disc = ldotc * ldotc - csq * lsq + lsq
            if disc < 0.0:
                disc = 0.0
            d1 = (ldotc - math.sqrt(disc)) / lsq
            ex, ey, ez = rx * d1 + px, ry * d1 + py, rz * d1 + pz

            lon_deg = math.atan2(ey, ex) * r2d
            p_xy = math.sqrt(ex * ex + ey * ey)
            za, pb = ez * a, p_xy * b
            r_t = math.sqrt(za * za + pb * pb)
            sin_t, cos_t = za / r_t, pb / r_t
            lat = math.atan2(ez + bow_ep2 * b * sin_t ** 3, p_xy - bow_e2 * a * cos_t ** 3)
            lat_out[i] = lat * r2d

            sin_lat, cos_lat = math.sin(lat), math.cos(lat)
            nrad = a / math.sqrt(1.0 - bow_e2 * sin_lat * sin_lat)
            if abs(sin_lat) > 0.9998476951563913:
                alt_out[i] = (abs(ez) / abs(sin_lat) - nrad * (1.0 - bow_e2)) * 1000.0
            else:
                alt_out[i] = (p_xy / cos_lat - nrad) * 1000.0

            lon_deg -= gmst[i] * r2d
            if lon_deg < -180.0:
                lon_deg += 360.0
            lon_out[i] = lon_deg

        return lon_out, lat_out, alt_out


def _interpolate_scan_endpoint_states(orb, times):
    """Interpolate short scan arcs from propagated endpoint states."""
    if isinstance(orb, (list, tuple)):
        orb = Orbital("mysatellite", line1=orb[0], line2=orb[1])
    endpoint_times = times[:, [0, -1]]
    endpoint_pos, endpoint_vel = orb.get_position(endpoint_times.reshape(-1), normalize=False)
    n_rows, n_pixels = times.shape
    positions = endpoint_pos.reshape(3, n_rows, 2)
    velocities = endpoint_vel.reshape(3, n_rows, 2)
    fraction = _scan_fraction(times, endpoint_times)
    pos = positions[:, :, :1] + np.diff(positions, axis=2) * fraction[np.newaxis]
    vel = velocities[:, :, :1] + np.diff(velocities, axis=2) * fraction[np.newaxis]
    return pos.reshape(3, n_rows * n_pixels), vel.reshape(3, n_rows * n_pixels)


def _scan_fraction(times, endpoints):
    """Return each pixel's fractional position between its scan endpoints."""
    duration = (endpoints[:, 1] - endpoints[:, 0]).astype("timedelta64[ns]").astype(float)
    elapsed = (times - endpoints[:, :1]).astype("timedelta64[ns]").astype(float)
    return np.divide(elapsed, duration[:, None], out=np.zeros_like(elapsed), where=duration[:, None] != 0)


def _geolocate_fused(orb, sgeom, times, roll, pitch, yaw=0.0, yaw_steering=False,
                     nadir_convention=None, rotation_order=None):
    """Inner fused-path implementation for :func:`geolocate`."""
    if isinstance(orb, (list, tuple)):
        tle1, tle2 = orb
        orb = Orbital("mysatellite", line1=tle1, line2=tle2)

    times_arr = np.asanyarray(times)
    pos, vel = _get_satpos(orb, times_arr, lines_per_scan=sgeom.lines_per_scan)

    M = pos.shape[1]          # number of scan rows
    N = sgeom.fovs.shape[2]   # pixels per row

    # Local orbital frame: uses geodetic nadir (subpoint) matching _local_frame
    nadir, _along, cross_track = _local_frame(pos, vel, nadir_convention)

    # Cross-track trig: same for every scan row → compute on N values only
    cos_c = np.cos(sgeom.fovs[0][0, :] + roll).astype(np.float64)    # (N,)
    sin_c = np.sin(sgeom.fovs[0][0, :] + roll).astype(np.float64)    # (N,)

    # Along-track angle: same for every pixel in a row → one value per row
    along_angles = (sgeom.fovs[1][:, 0] + pitch).astype(np.float64)  # (M,)
    if np.any(along_angles):
        rotation_order = _resolve_convention(
            rotation_order, "rotation_order", "rotation order", "rotation_order='pitch_first'"
        )
    pitch_first = rotation_order != "legacy"

    # Effective yaw per scan line: explicit bias + optional steering term
    if yaw_steering:
        steering = compute_yaw_steering(pos, vel)          # (M,) array
        yaw_angles = (yaw + steering).astype(np.float64)
    else:
        yaw_angles = np.full(M, yaw, dtype=np.float64)

    # GMST once per scan row
    if times_arr.ndim == 2:
        gmst_scan = astronomy.gmst(times_arr[:, 0]).astype(np.float64)
    else:
        gmst_scan = astronomy.gmst(times_arr[::N][:M]).astype(np.float64)

    lon_out = np.empty(M * N)
    lat_out = np.empty(M * N)
    alt_out = np.empty(M * N)

    _fused_geolocate_numba(
        np.ascontiguousarray(pos, dtype=np.float64),
        np.ascontiguousarray(nadir, dtype=np.float64),
        np.ascontiguousarray(cross_track, dtype=np.float64),
        cos_c, sin_c, along_angles, pitch_first, yaw_angles, gmst_scan,
        lon_out, lat_out, alt_out,
    )
    return lon_out, lat_out, alt_out


def _gmst_per_pixel(times):
    """Return per-pixel GMST in radians, computing once per scan line for 2-D arrays.

    For 2-D time arrays (shape ``(N_scans, N_pixels)``) GMST is evaluated only
    for the first pixel time of each scan line, then broadcast back to the full
    pixel count.  This reduces the computation from ``N_scans * N_pixels``
    evaluations down to ``N_scans``.
    """
    times = np.asanyarray(times)
    if times.ndim == 2:
        return np.repeat(astronomy.gmst(times[:, 0]), times.shape[1])
    return astronomy.gmst(times)


def _get_satpos(orb, times, lines_per_scan=1):
    """Compute satellite position/velocity, calling SGP4 once per scan line.

    **Context**: the original pyorbital used a pure-numpy SGP4 propagator
    that accepted arrays of any shape, so passing a 2-D ``(N_scans, N_pixels)``
    time array worked transparently and produced per-pixel positions.  The
    current implementation uses the ``sgp4`` C library whose ``sgp4_array``
    only accepts 1-D input, so 2-D time arrays must be handled explicitly.

    **Approximation**: for 2-D time arrays (shape ``(N_scans, N_pixels)``),
    SGP4 is evaluated only for the first pixel time of each scan.  The
    satellite moves at roughly 7.5 km/s; across a typical instrument scan
    (< 0.5 s) this amounts to a position error of < 4 km, producing a
    geolocation error well below 0.01° — smaller than the accuracy of the
    sample-and-interpolate approach it replaces.

    For multi-line instruments (``lines_per_scan > 1``), the rows time array
    has shape ``(N_scans * lines_per_scan, N_pixels)``.  SGP4 is evaluated
    once per *instrument* scan (striding by ``lines_per_scan``), and the
    resulting per-scan positions are repeated to cover all output rows.

    For 1-D time arrays the usual per-element SGP4 is used unchanged.
    """
    times = np.asanyarray(times)
    if times.ndim == 2:
        scan_times = times[::lines_per_scan, 0]
        pos, vel = orb.get_position(scan_times, normalize=False)
        if lines_per_scan > 1:
            pos = np.repeat(pos, lines_per_scan, axis=1)
            vel = np.repeat(vel, lines_per_scan, axis=1)
        return pos, vel
    return orb.get_position(times, normalize=False)


def compute_pixels(orb, sgeom, times, rpy=(0.0, 0.0, 0.0), yaw_steering=False,
                   nadir_convention=None, rotation_order=None):
    """Compute cartesian coordinates of the pixels in instrument scan."""
    if isinstance(orb, (list, tuple)):
        tle1, tle2 = orb
        orb = Orbital("mysatellite", line1=tle1, line2=tle2)

    pos, vel = _get_satpos(orb, times, lines_per_scan=sgeom.lines_per_scan)

    vectors = sgeom.vectors(pos, vel, *rpy, yaw_steering=yaw_steering,
                            nadir_convention=nadir_convention,
                            rotation_order=rotation_order)

    # Compute intersection of pixel lines with the ellipsoid declared by A and B.
    # http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    radius = _ellipsoid_radius_scaling()

    n_total = vectors.shape[1]
    n_pos = pos.shape[1]

    if n_pos == n_total:
        return _ellipsoid_intersection(vectors, pos, radius)

    # pos has fewer columns than vectors (one per scan line, many pixels each).
    # Avoid np.repeat which would allocate a full (3, n_total) copy of pos.
    pixels_per_line = n_total // n_pos
    return _ellipsoid_intersection_broadcast(vectors, pos, radius, pixels_per_line)



def _ellipsoid_intersection(vectors, pos, radius):
    """Intersect pixel vectors with the WGS-84 ellipsoid (pos already broadcast).

    When a ray misses the ellipsoid (discriminant < 0) the discriminant is
    clamped to 0, returning the point of closest approach on the tangent line.
    This keeps results finite and continuous during optimizer exploration.
    """
    centre = -pos
    shape = vectors.shape
    xr_ = vectors.reshape([3, -1]) * radius
    cr_ = centre.reshape([3, -1]) * radius
    ldotc = np.einsum("ij,ij->j", xr_, cr_)
    lsq = np.einsum("ij,ij->j", xr_, xr_)
    csq = np.einsum("ij,ij->j", cr_, cr_)
    discriminant = np.maximum(ldotc ** 2 - csq * lsq + lsq, 0.0)
    d1_ = (ldotc - np.sqrt(discriminant)) / lsq
    return vectors * d1_.reshape(shape[1:]) - centre.reshape(shape)


def _ellipsoid_intersection_broadcast(vectors, pos, radius, pixels_per_line):
    """Intersect pixel vectors with WGS-84 ellipsoid using per-scan-line broadcasting.

    Avoids allocating a full (3, n_total) copy of *pos* when many pixels share
    the same satellite position (i.e. when pos has one column per scan line and
    vectors has *pixels_per_line* columns per scan line).

    Args:
        vectors: Pixel unit vectors, shape ``(3, n_lines * pixels_per_line)``.
        pos: Satellite positions, shape ``(3, n_lines)``.
        radius: WGS-84 axis scaling vector, shape ``(3, 1)``.
        pixels_per_line: Number of pixels per scan line.

    Returns:
        Cartesian surface coordinates, shape ``(3, n_lines * pixels_per_line)``.
    """
    n_lines = pos.shape[1]
    # Reshape vectors to (3, n_lines, pixels_per_line) — no copy, just a view.
    vec3 = vectors.reshape(3, n_lines, pixels_per_line)
    xr_ = vec3 * radius.reshape(3, 1, 1)                    # (3, M, K)
    cr_ = -pos * radius                                       # (3, M)
    ldotc = np.einsum("ijk,ij->jk", xr_, cr_)               # (M, K)
    lsq = np.einsum("ijk,ijk->jk", xr_, xr_)                # (M, K)
    csq = np.einsum("ij,ij->j", cr_, cr_)                   # (M,)
    csq_exp = csq[:, np.newaxis]
    discriminant = np.maximum(ldotc ** 2 - csq_exp * lsq + lsq, 0.0)
    d1_ = (ldotc - np.sqrt(discriminant)) / lsq
    result = vec3 * d1_[np.newaxis] + pos[:, :, np.newaxis]
    return result.reshape(3, -1)


def norm(v):
    """Return the norm of the vector *v*."""
    return np.sqrt(np.dot(v, v.conj()))


def mnorm(m, axis=None):
    """Norm of a matrix of vectors stacked along the *axis* dimension."""
    if axis is None:
        axis = np.ndim(m) - 1
    return np.sqrt((m**2).sum(axis))


def vnorm(m):
    """Norms of a matrix of column vectors."""
    return np.sqrt((m**2).sum(0))


def hnorm(m):
    """Norms of a matrix of row vectors."""
    return np.sqrt((m**2).sum(1))


def _n_scans_from_duration(duration, time_sampling):
    """Compute the number of scan lines covering a given duration."""
    if time_sampling <= np.timedelta64(0, "ns"):
        return 1
    return max(1, int(duration / time_sampling))


def _sample_indices(count, n_samples):
    """Return n_samples evenly-spaced integer indices from 0 to count-1."""
    return np.linspace(0, count - 1, n_samples, dtype=int)


def _boundary_scan_pixel_pairs(scan_indices, pixel_indices):
    """Build ordered (scan, pixel) pairs tracing the swath boundary polygon.

    Traverses: bottom edge left→right, right edge bottom→top,
    top edge right→left, left edge top→bottom, closing point.
    """
    scans, pixels = [], []

    for px in pixel_indices:               # bottom edge: left → right
        scans.append(scan_indices[0])
        pixels.append(px)
    for sc in scan_indices[1:]:            # right edge: bottom → top
        scans.append(sc)
        pixels.append(pixel_indices[-1])
    for px in reversed(pixel_indices[:-1]): # top edge: right → left
        scans.append(scan_indices[-1])
        pixels.append(px)
    for sc in reversed(scan_indices[1:-1]): # left edge: top → bottom
        scans.append(sc)
        pixels.append(pixel_indices[0])

    scans.append(scans[0])     # close the polygon
    pixels.append(pixels[0])
    return np.array(scans), np.array(pixels)


def _geolocate_boundary(swath, edge_scans, edge_pixels, start_time, tle, rpy):
    """Geolocate a set of (scan, pixel) boundary pairs and return (lons, lats)."""
    x_fovs, y_fovs = swath.scanline.angles(edge_pixels)
    sgeom = ScanGeometry(np.vstack((x_fovs, y_fovs)), edge_scans * swath.time_sampling)
    s_times = sgeom.times(start_time)
    pixels_pos = compute_pixels(tle, sgeom, s_times, rpy)
    lons, lats, _ = get_lonlatalt(pixels_pos, s_times)
    return lons, lats


def bounding_box(swath, start_time, end_time, tle, points_per_edge=10, rpy=(0.0, 0.0, 0.0)):
    """Compute a bounding polygon for a satellite swath.

    Args:
        swath: A PushbroomSwath with scanline and time_sampling attributes.
        start_time: Start of the observation period (datetime).
        end_time: End of the observation period (datetime).
        tle: TLE data as (line1, line2) tuple.
        points_per_edge: Number of sample points per edge (including corners).
        rpy: Roll, pitch, yaw corrections.

    Returns:
        Tuple of (lons, lats) arrays forming a closed polygon.
    """
    duration = np.datetime64(end_time) - np.datetime64(start_time)
    n_scans = _n_scans_from_duration(duration, swath.time_sampling)

    scan_indices = _sample_indices(n_scans, points_per_edge)
    pixel_indices = _sample_indices(swath.scanline.pixels_per_scan, points_per_edge)

    edge_scans, edge_pixels = _boundary_scan_pixel_pairs(scan_indices, pixel_indices)
    return _geolocate_boundary(swath, edge_scans, edge_pixels, start_time, tle, rpy)


if __name__ == "__main__":
    # NOAA 18 (from the 2011-10-12, 16:55 utc)
    # 1 28654U 05018A   11284.35271227  .00000478  00000-0  28778-3 0  9246
    # 2 28654  99.0096 235.8581 0014859 135.4286 224.8087 14.11526826329313

    noaa18_tle1 = "1 28654U 05018A   11284.35271227  .00000478  00000-0  28778-3 0  9246"
    noaa18_tle2 = "2 28654  99.0096 235.8581 0014859 135.4286 224.8087 14.11526826329313"

    import datetime as dt
    t = dt.datetime(2011, 10, 12, 13, 45)

    # edge and centre of an avhrr scanline
    # sgeom = ScanGeometry([(-0.9664123687741623, 0),
    #                      (0, 0)],
    #                     [0, 0.0, ])
    # print compute_pixels((noaa18_tle1, noaa18_tle2), sgeom, t)

    # avhrr swath
    scanline_nb = 1

    # building the avhrr angles, 2048 pixels from +55.37 to -55.37 degrees
    avhrr = np.vstack(((np.arange(2048) - 1023.5) / 1024 * np.deg2rad(-55.37),
                       np.zeros((2048,)))).transpose()
    avhrr = np.tile(avhrr, [scanline_nb, 1])
    # building the corresponding times array
    offset = np.arange(scanline_nb) * 0.1667
    times = (np.tile(np.arange(2048) * 0.000025 + 0.0025415, [scanline_nb, 1])
             + np.expand_dims(offset, 1))
    # build the scan geometry object
    sgeom = ScanGeometry(avhrr, times.ravel())

    # print the lonlats for the pixel positions
    s_times = sgeom.times(t)
    pixels_pos = compute_pixels((noaa18_tle1, noaa18_tle2), sgeom, s_times)
    print(get_lonlatalt(pixels_pos, s_times))
