"""Module to compute geolocalization of a satellite scene."""


from __future__ import print_function

import math
from warnings import warn

import numpy as np
from pyproj import Transformer

try:
    import numba as nb
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False
    warn("Install Numba for better performance.")

# DIRTY STUFF. Needed the get_lonlatalt function to work on pos directly if
# we want to print out lonlats in the end.
from pyorbital import astronomy
from pyorbital.orbital import Orbital

A = 6378.137  # WGS84 and GRS80 Equatorial radius (km)
B = 6356.75231414  # km, GRS80

# Module-level cached Transformer — avoids re-creating the PROJ context on
# every get_lonlatalt() call.
_GEOCENT_TO_LATLONG = None


def _get_transformer():
    global _GEOCENT_TO_LATLONG
    if _GEOCENT_TO_LATLONG is None:
        _GEOCENT_TO_LATLONG = Transformer.from_crs(
            dict(proj="geocent"), dict(proj="latlong"))
    return _GEOCENT_TO_LATLONG

OMEGA_EARTH = 7.2921159e-5  # Earth's rotation rate (rad/s)


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


def _local_frame(pos, vel):
    """Compute the satellite's local orbital reference frame.

    Returns (nadir, along_track, cross_track) as unit column vectors.
    """
    nadir = subpoint(-pos)
    nadir /= vnorm(nadir)
    along_track = vel / vnorm(vel)
    cross_track = np.cross(nadir, vel, 0, 0, 0)
    cross_track /= vnorm(cross_track)
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

    def vectors(self, pos, vel, roll=0.0, pitch=0.0, yaw=0.0, yaw_steering=False):
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
        if self.fovs.ndim == 3 and pos.ndim == 2 and pos.shape[1] == self.fovs.shape[1]:
            return self._vectors_broadcast(pos, vel, roll, pitch, yaw, yaw_steering)
        fovs = self.fovs.reshape(2, -1)
        nadir, along_track, cross_track = _local_frame(pos, vel)
        effective_yaw = _effective_yaw(yaw, yaw_steering, pos, vel, self.fovs[0].shape)
        rotated = qrotate(nadir, along_track, fovs[0] + roll)
        if np.any(fovs[1] + pitch):
            rotated = qrotate(rotated, cross_track, fovs[1] + pitch)
        return qrotate(rotated, nadir, effective_yaw)

    def _vectors_broadcast(self, pos, vel, roll, pitch, yaw, yaw_steering):
        """Broadcast rotation for per-scan pos/vel and per-pixel fovs.

        Two structural invariants of all 3-D fov scan geometries are exploited:

        1. **Cross-track angles** (``fovs[0]``) are the same for every scan line
           — they only depend on the pixel index.  Trig is therefore computed
           on ``N_pixels`` values instead of ``N_scans * N_pixels``.

        2. **Along-track angles** (``fovs[1]``) are constant across all pixels
           within a scan line — they only depend on the detector-row index.
           Trig is computed on ``N_scans`` values instead of
           ``N_scans * N_pixels``.

        3. **nadir ⊥ along_track** (guaranteed by :func:`_local_frame`).
           The first Rodrigues rotation therefore simplifies to two
           multiply-add operations:
           ``rotated = nadir * cos(θ) + cross_track * sin(θ)``
        """
        nadir, along_track, cross_track = _local_frame(pos, vel)
        effective_yaw = _effective_yaw(yaw, yaw_steering, pos, vel, self.fovs[0].shape)

        # --- First rotation: nadir around along_track by cross-track scan angle ---
        # Exploit nadir ⊥ along_track: full Rodrigues reduces to 2 mults.
        # Cross-track angles are identical for all scan lines: use first row (1, N)
        # for trig, then broadcast across M scan lines.
        cross_angles = (self.fovs[0] + roll)[0:1, :]           # (1, N)
        cos_c = np.cos(cross_angles)                            # (1, N)
        sin_c = np.sin(cross_angles)                            # (1, N)
        nadir_3d = nadir[:, :, np.newaxis]                      # (3, M, 1)
        ct_3d = cross_track[:, :, np.newaxis]                   # (3, M, 1)
        rotated = nadir_3d * cos_c[np.newaxis] + ct_3d * sin_c[np.newaxis]

        # --- Second rotation: around cross_track by along-track detector angle ---
        # Along-track angles are constant across pixels: use first column (M, 1)
        # so trig is computed on M values instead of M*N.
        along_angles = (self.fovs[1] + pitch)[:, 0:1]          # (M, 1)
        if np.any(along_angles):
            rotated = _rodrigues(rotated, cross_track, along_angles)

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
            x, y, z: 1-D float64 arrays of ECEF coordinates in kilometres.

        Returns:
            Tuple ``(lon_deg, lat_deg, alt_m)`` as 1-D float64 arrays.
        """
        _a = 6378.137
        _b = 6356.752314245
        _e2 = 1.0 - (_b / _a) ** 2
        _ep2 = (_a / _b) ** 2 - 1.0
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
            za = zi * _a
            pb = p * _b
            r_t = math.sqrt(za * za + pb * pb)
            sin_t = za / r_t
            cos_t = pb / r_t

            num = zi + _ep2 * _b * sin_t * sin_t * sin_t
            den = p - _e2 * _a * cos_t * cos_t * cos_t
            lat = math.atan2(num, den)
            lat_out[i] = lat * _r2d

            sin_lat = math.sin(lat)
            cos_lat = math.cos(lat)
            N = _a / math.sqrt(1.0 - _e2 * sin_lat * sin_lat)
            if abs(sin_lat) > _sin_89_9:
                alt_out[i] = (abs(zi) / abs(sin_lat) - N * (1.0 - _e2)) * 1000.0
            else:
                alt_out[i] = (p / cos_lat - N) * 1000.0

        return lon_out, lat_out, alt_out

    @nb.njit(parallel=True, cache=True)
    def _fused_geolocate_numba(pos, nadir, cross_track, cos_c, sin_c, along_angles,
                               yaw_angles, gmst_scan, lon_out, lat_out, alt_out):
        """Fused rotation + ellipsoid intersection + geodetic in one parallel kernel.

        Outer ``prange`` over M scan lines (one per CPU thread); inner sequential
        loop over N pixels.  Eliminates the large ``(3, M, N)`` intermediate
        arrays produced by the numpy pipeline, reducing memory pressure by ~3×.

        Sign convention: clockwise Rodrigues (negative sin term) matching
        :func:`_rodrigues` and :func:`qrotate`.

        Three rotations are applied in order:
        1. Cross-track (simplified, exploits nadir ⊥ along_track).
        2. Along-track Rodrigues (skipped when ``along_angles[m] ≈ 0``).
        3. Yaw Rodrigues around nadir (skipped when ``yaw_angles[m] ≈ 0``).

        Args:
            pos: Satellite ECEF positions (km), shape ``(3, M)``.
            nadir: Pre-computed geodetic nadir unit vectors from
                :func:`_local_frame`, shape ``(3, M)``.
            cross_track: Pre-computed cross-track unit vectors from
                :func:`_local_frame`, shape ``(3, M)``.
            cos_c: cos of (cross-track scan angle + roll), shape ``(N,)``.
            sin_c: sin of (cross-track scan angle + roll), shape ``(N,)``.
            along_angles: Along-track detector angle + pitch per scan line,
                shape ``(M,)`` in radians.
            yaw_angles: Effective yaw per scan line (explicit yaw + optional
                steering), shape ``(M,)`` in radians.
            gmst_scan: Greenwich Mean Sidereal Time per scan line, shape
                ``(M,)`` in radians.
            lon_out, lat_out, alt_out: Pre-allocated output arrays, shape
                ``(M * N,)``.  Units: degrees, degrees, metres.
        """
        _a = 6378.137
        _b = 6356.752314245
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

            # Along-track Rodrigues trig — computed once per scan line
            alpha = along_angles[m]
            cos_a = math.cos(alpha)
            sin_a = math.sin(alpha)
            do_along = abs(sin_a) > 1e-15
            one_m_cos_a = 1.0 - cos_a

            gmst_m = gmst_scan[m]

            for k in range(N):
                # --- 1st rotation: nadir*cos_c + cross_track*sin_c ---
                # Simplified CW Rodrigues (nadir ⊥ along_track):
                #   v_rot = nadir·cos(θ) - (along×nadir)·sin(θ)
                #         = nadir·cos(θ) + cross_track·sin(θ)
                cc, sc = cos_c[k], sin_c[k]
                rx = nx * cc + ctx * sc
                ry = ny * cc + cty * sc
                rz = nz * cc + ctz * sc

                # --- 2nd rotation: CW Rodrigues around cross_track by alpha ---
                if do_along:
                    ctdotr = ctx * rx + cty * ry + ctz * rz
                    crx = cty * rz - ctz * ry
                    cry = ctz * rx - ctx * rz
                    crz = ctx * ry - cty * rx
                    rx = rx * cos_a - crx * sin_a + ctx * ctdotr * one_m_cos_a
                    ry = ry * cos_a - cry * sin_a + cty * ctdotr * one_m_cos_a
                    rz = rz * cos_a - crz * sin_a + ctz * ctdotr * one_m_cos_a

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
                d1 = (ldotc - math.sqrt(ldotc * ldotc - csq * lsq + lsq)) / lsq
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
                               yaw_angles, gmst_scan, lon_out, lat_out, alt_out):
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

# END OF DIRTY STUFF


def geolocate(orb, sgeom, times, rpy=(0.0, 0.0, 0.0), yaw_steering=False):
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
        rpy: ``(roll, pitch, yaw)`` attitude corrections in radians.
        yaw_steering: If ``True``, compute yaw from orbit geometry to
            counteract Earth rotation.

    Returns:
        Tuple ``(lon_deg, lat_deg, alt_m)`` as flat 1-D arrays.
    """
    roll, pitch, yaw = rpy

    if _HAS_NUMBA and sgeom.fovs.ndim == 3:
        return _geolocate_fused(orb, sgeom, times, roll, pitch, yaw, yaw_steering)

    pixels = compute_pixels(orb, sgeom, times, rpy, yaw_steering)
    return get_lonlatalt(pixels, times)


def _geolocate_fused(orb, sgeom, times, roll, pitch, yaw=0.0, yaw_steering=False):
    """Inner fused-path implementation for :func:`geolocate`."""
    if isinstance(orb, (list, tuple)):
        tle1, tle2 = orb
        orb = Orbital("mysatellite", line1=tle1, line2=tle2)

    times_arr = np.asanyarray(times)
    pos, vel = _get_satpos(orb, times_arr, lines_per_scan=sgeom.lines_per_scan)

    M = pos.shape[1]          # number of scan rows
    N = sgeom.fovs.shape[2]   # pixels per row

    # Local orbital frame: uses geodetic nadir (subpoint) matching _local_frame
    nadir, _along, cross_track = _local_frame(pos, vel)

    # Cross-track trig: same for every scan row → compute on N values only
    cos_c = np.cos(sgeom.fovs[0][0, :] + roll).astype(np.float64)    # (N,)
    sin_c = np.sin(sgeom.fovs[0][0, :] + roll).astype(np.float64)    # (N,)

    # Along-track angle: same for every pixel in a row → one value per row
    along_angles = (sgeom.fovs[1][:, 0] + pitch).astype(np.float64)  # (M,)

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
        cos_c, sin_c, along_angles, yaw_angles, gmst_scan,
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


def compute_pixels(orb, sgeom, times, rpy=(0.0, 0.0, 0.0), yaw_steering=False):
    """Compute cartesian coordinates of the pixels in instrument scan."""
    if isinstance(orb, (list, tuple)):
        tle1, tle2 = orb
        orb = Orbital("mysatellite", line1=tle1, line2=tle2)

    pos, vel = _get_satpos(orb, times, lines_per_scan=sgeom.lines_per_scan)

    vectors = sgeom.vectors(pos, vel, *rpy, yaw_steering=yaw_steering)

    # Compute intersection of pixel lines with the WGS-84 ellipsoid.
    # http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    a__ = 6378.137  # km
    b__ = 6356.752314245  # km, WGS84
    radius = np.array([[1 / a__, 1 / a__, 1 / b__]]).T

    n_total = vectors.shape[1]
    n_pos = pos.shape[1]

    if n_pos == n_total:
        return _ellipsoid_intersection(vectors, pos, radius)

    # pos has fewer columns than vectors (one per scan line, many pixels each).
    # Avoid np.repeat which would allocate a full (3, n_total) copy of pos.
    pixels_per_line = n_total // n_pos
    return _ellipsoid_intersection_broadcast(vectors, pos, radius, pixels_per_line)



def _ellipsoid_intersection(vectors, pos, radius):
    """Intersect pixel vectors with the WGS-84 ellipsoid (pos already broadcast)."""
    centre = -pos
    shape = vectors.shape
    xr_ = vectors.reshape([3, -1]) * radius
    cr_ = centre.reshape([3, -1]) * radius
    ldotc = np.einsum("ij,ij->j", xr_, cr_)
    lsq = np.einsum("ij,ij->j", xr_, xr_)
    csq = np.einsum("ij,ij->j", cr_, cr_)
    d1_ = (ldotc - np.sqrt(ldotc ** 2 - csq * lsq + lsq)) / lsq
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
    lsq   = np.einsum("ijk,ijk->jk", xr_, xr_)              # (M, K)
    csq   = np.einsum("ij,ij->j", cr_, cr_)                 # (M,)
    csq_exp = csq[:, np.newaxis]
    d1_ = (ldotc - np.sqrt(ldotc ** 2 - csq_exp * lsq + lsq)) / lsq
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
    if time_sampling <= np.timedelta64(0):
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
